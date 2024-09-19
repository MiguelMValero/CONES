from IPython.core import magic_arguments
from IPython.core.magic import cell_magic, line_magic, Magics, magics_class
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
import numpy as np 

def mag2(a):
  return np.dot(a, a)

def distance2_point_point(a, b):
  return mag2(a-b)

def distance2_point_line(p, a, b):
    v = b - a
    t = np.dot(p - a, v) / mag2(v)
    if t < 0:
        c = a
    elif t > 1:
        c = b
    else:
        c = (1-t)*a + t*b
    return c, mag2(p - c)


def ray_tracing(x,y,poly):
    n = len(poly)
    inside = False
    p2x = 0.0
    p2y = 0.0
    xints = 0.0
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def point_location(tolerance, method):
    fig, ax = plt.subplots()

    # polygon
    poly_coord = np.array([
        [3.36, 1.85],
        [-2.1, 2.75],
        [-4.88, -1.13],
        [-0.34, -4.21],
        [4.14, -3.03],
        [1.82, -0.79]])
    n_vtx = len(poly_coord)

    # points
    pts_coord = np.array([
        [-0.3, -1.6],
        [0.68, 3.6],
        [-3.7, -3.3]])
    n_pts = len(pts_coord)

    scale = 0.75
    poly_coord *= scale
    pts_coord  *= scale

    # bbox
    xymin = np.amin(poly_coord, axis=0)
    xymax = np.amax(poly_coord, axis=0)

    xyeps = tolerance*(xymax - xymin)
    xymin -= xyeps
    xymax += xyeps


    if method in ["OCTREE", "DBBTREE"]:
        is_located = np.zeros(n_pts, dtype=bool)
        for i, p in enumerate(pts_coord):
            is_located[i] = p[0] >= xymin[0] and p[0] <= xymax[0] and p[1] >= xymin[1] and p[1] <= xymax[1]
    else:
        is_located = np.ones(n_pts, dtype=bool)

    proj_coord = np.empty(pts_coord.shape, dtype=float)
    for i, p in enumerate(pts_coord):
        if is_located[i]:
            if ray_tracing(p[0], p[1], poly_coord):
                proj_coord[i] = p
            else:
                dmin = 1e30
                for j in range(n_vtx):
                    k = (j+1)%n_vtx
                    cp, d = distance2_point_line(p, poly_coord[j], poly_coord[k])
                    if d < dmin:
                        dmin = d
                        proj_coord[i] = cp

    # plot
    ax.add_patch(Polygon(poly_coord,
                         fill=True,
                         fc="#f0f8ff",
                         ec="k"))

    if method in ["OCTREE", "DBBTREE"]:
        ax.add_patch(Rectangle(xymin, xymax[0] - xymin[0], xymax[1] - xymin[1],
                               fill=False,
                               ec="0.7",
                               ls="--"))


    for i, p in enumerate(pts_coord):
        if is_located[i]:
            ax.plot([p[0], proj_coord[i,0]], [p[1], proj_coord[i,1]], "g:")
            ax.plot(proj_coord[i,0], proj_coord[i,1], "go", ms=7, mfc="w")

        ax.plot(p[0], p[1], "o", color="g" if is_located[i] else "r", ms=4)

    # fake pts for legend
    ax.plot(100, 100, "ro", ms=4,          label="unlocated")
    ax.plot(100, 100, "go", ms=4,          label="located")
    ax.plot(100, 100, "go", ms=7, mfc="w", label="projection")

    ax.legend(loc='upper left', edgecolor="w")

    ax.set_aspect("equal")
    ax.set_xlim(-6,6)
    ax.set_ylim(-6,4)
    plt.axis("off")
    plt.show()

@magics_class
class FigureMagics(Magics):
    """
    Plot interactive figures
    """
    @line_magic
    @magic_arguments.magic_arguments()

    def localization(self, line):
        """
        Localization figure
        """

        from ipywidgets import interactive, FloatSlider, Dropdown, HBox
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle, Polygon
        import numpy as np

        tol = FloatSlider(min=0., max=0.3, step=1e-4, value=0, continuous_update=True)
        met = Dropdown(options=["OCTREE", "LOCATE_ALL_TGT"], value="OCTREE")
        ui = HBox([tol, met])

        interactive_plot = interactive(point_location, tolerance=tol, method=met)
        out = interactive_plot.children[-1]
        out.layout.height = '300px'
        display(ui, out)

def load_ipython_extension(ipython):
  ipython.register_magics(FigureMagics)
