from IPython.core import magic_arguments
from IPython.core.magic import cell_magic, line_magic, Magics, magics_class

def nice_view(plotter):
  import numpy as np
  p_bounds = np.asarray(plotter.bounds)
  p_range = p_bounds[1::2] - p_bounds[0::2]
  imin = np.argmin(p_range)

  if p_range[imin] < 1e-14:
    if imin == 0:
      plotter.view_yz()
    elif imin == 1:
      plotter.view_zx()
    else:
      plotter.view_xy()

def read_geom(filename):
  import re
  import numpy as np
  coord    = None
  face_vtx = []
  n_elt_type = dict()

  with open(filename, "r") as f:
    while True:
      line = f.readline()
      if line == "": break

      match = re.search("coordinates", line)
      if match is not None:
        n_vtx = int(f.readline())
        points = np.empty((n_vtx, 3), dtype=float)
        for j in range(3):
          for i in range(n_vtx):
            points[i,j] = float(f.readline())
        coord = points

      for kw in ["point", "tria3", "quad4"]:
        match = re.search(kw, line)
        if match is not None:
          n_elt = int(f.readline())
          n_elt_type[kw] = n_elt
          for i in range(n_elt):
            vertices = [int(a)-1 for a in f.readline().split()]
            face_vtx.append(vertices)
          break

      match = re.search("nsided", line)
      if match is not None:
        n_elt = int(f.readline())
        n_elt_type["nsided"] = n_elt
        for i in range(n_elt):
          f.readline()

        for i in range(n_elt):
          vertices = [int(a)-1 for a in f.readline().split()]
          face_vtx.append(vertices)

  return {"points": coord, "face_vtx": face_vtx, "n_elt_type": n_elt_type}


def read_node_scalar(filename, n_node):
  import re
  import numpy as np
  values = None
  with open(filename, "r") as f:
    while True:
      line = f.readline()
      if line == "": break
      match = re.search("coordinates", line)
      if match is not None:
        values = np.empty(n_node, dtype=float)
        for i in range(n_node):
          values[i] = float(f.readline())

  return values

def read_element_scalar(filename, n_elt_type):
  import re
  import numpy as np
  tn_elt = 0
  for kw in n_elt_type:
    tn_elt += n_elt_type[kw]
  values = np.empty(tn_elt, dtype=float)
  idx = 0
  with open(filename, "r") as f:
    while True:
      line = f.readline()
      if line == "": break

      for kw in n_elt_type:
        match = re.search(kw, line)
        if match is not None:
          n_elt = n_elt_type[kw]
          for i in range(n_elt):
            values[idx] = float(f.readline())
            idx += 1

  return values




def visu_n_files(files,
                 fields=[""],
                 styles=[""],
                 same_view=False,
                 lighting=True,
                 show_grid=False):
  try:
    import pyvista as pv
    import numpy   as np

    n_files = len(files)
    if n_files == 1:
      same_view = True

    # pyvista options
    pv.set_plot_theme("document")
    pv.global_theme.trame.interactive_ratio = 2
    pv.global_theme.trame.still_ratio       = 2


    # plotter
    if same_view:
      p = pv.Plotter(notebook=True)
    else:
      n_col = 2 # 2 subplots per column
      n_row = n_files//n_col
      p = pv.Plotter(notebook=True, shape=(n_row, n_col))

    p.background_color = 'w'
    p.enable_anti_aliasing()

    for i_file, filename in enumerate(files):
      if not same_view:
        i_col = i_file%n_col
        i_row = i_file//n_col
        p.subplot(i_row, i_col)

      # load mesh file
      try:
        mesh = pv.read(filename)
        if not isinstance(mesh, pv.MultiBlock):
          mesh = [mesh]

        for block in mesh:

          scalars = None
          if i_file < len(fields):
            if len(fields[i_file]) > 0:
              scalars = fields[i_file]

          style = None
          if i_file < len(styles):
            if len(styles[i_file]) > 0:
              style = styles[i_file]

          if scalars is None:
            color = [0.8]*3
          else:
            color = None

          p.add_mesh(block,
                     show_edges=True,
                     style=style,
                     cmap="coolwarm",
                     color=color,
                     scalars=scalars,
                     lighting=lighting,
                     point_size=6)
      except:
        print(f"Output file '{filename}' not found", flush=True)

    p.link_views()

    nice_view(p)

    # show loaded meshes
    p.show(jupyter_backend='client')

  except ImportError:
    print("module not found: pyvista (required for interactive visualization)")


@magics_class
class VisuMagics(Magics):
  """
  Visualize multiple files with pyvista
  """
  @cell_magic
  @magic_arguments.magic_arguments()
  @magic_arguments.argument('--same_view', '-sv',
                            help='Show all files in same view',
                            action="store_true")
  @magic_arguments.argument('--no_lighting', '-nl',
                            help='No lighting',
                            action="store_true")
  def visualize(self, line='', cell=None):
    args = magic_arguments.parse_argstring(self.visualize, line)

    # parse cell content
    cell_lines = cell.split("\n")
    n = len(cell_lines)
    files  = []
    fields = []
    styles = []

    for l in cell_lines:
      a = [k.rstrip().lstrip() for k in l.split(":")]

      if len(a[0]) > 0:
        files.append(a[0])

        if len(a) > 1:
          fields.append(a[1])
        else:
          fields.append("")

        if len(a) > 2:
          styles.append(a[2])
        else:
          styles.append("")

    visu_n_files(files     = files,
                 fields    = fields,
                 styles    = styles,
                 same_view = args.same_view,
                 lighting  = not args.no_lighting)


def load_ipython_extension(ipython):
  ipython.register_magics(VisuMagics)
