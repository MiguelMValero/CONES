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


def read_ensight(filename):
  import os
  import re
  import numpy as np
  import pyvista as pv
  from pyvista import CellType

  match_number = re.compile('\ *-?[0-9]+\.?[0-9]*(?:[Ee]\ *-?\+?\ *[0-9]+)?')

  path = os.path.dirname(filename)
  if len(path) > 0:
    if path[-1] != "/":
      path += "/"

  # read case
  geom_prefix = None
  n_steps = -1
  start   = -1
  incr    = -1
  with open(filename, "r") as f:
    while True:
      line = f.readline()
      if line == "": break

      if re.search("GEOMETRY", line) is not None:
        while True:
          line = f.readline()
          if line == "" or line == "\n": break
          match = re.search("model:", line)
          if match is not None:
            info = line[match.span()[1]:].split()
            idx = info[1].rfind(".")
            geom_prefix = info[1][:idx]
            geom_padding = len(info[1]) - idx - 1
            #geom_mode = info[2]
            break

      if re.search("VARIABLE", line) is not None:
        scalars = dict()
        while True:
          line = f.readline()
          if line == "" or line == "\n": break
          match = re.search("scalar per node:", line)
          if match is not None:
            info = line[match.span()[1]:].split()
            name = info[1]
            idx = info[2].rfind(".")
            prefix = info[2][:idx]
            padding = len(info[2]) - idx - 1
            scalars[name] = {"prefix": prefix, "padding": padding}

      if re.search("TIME", line) is not None:
        while True:
          line = f.readline()
          if line == "": break

          match = re.search("number of steps:", line)
          if match is not None:
            n_steps = int(line[match.span()[1]:])
            #print(f"n_steps = {n_steps}")

          match = re.search("filename start number:", line)
          if match is not None:
            start = int(line[match.span()[1]:])
            #print(f"start = {start}")

          match = re.search("filename increment:", line)
          if match is not None:
            incr = int(line[match.span()[1]:])
            #print(f"incr = {incr}")
            break

          # read time values?

  # read geometry
  instants = dict()
  #geom_prefix = "chr.geo"
  for i_step in range(n_steps):
    fmt = "{:0%d}" % geom_padding
    geom_name = path + geom_prefix + "." + fmt.format(start + incr*i_step)

    grid = None
    with open(geom_name, "r") as f:
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

        match = re.search("nsided", line)
        if match is not None:
          n_elt = int(f.readline())
          nside = np.empty(n_elt, dtype=int)
          size = 0
          for i in range(n_elt):
            nside[i] = int(f.readline())
            size += nside[i] + 1

          connec = np.empty(size, dtype=int)
          face_vtx = []
          idx = 0
          for i in range(n_elt):
            connec[idx] = nside[i]
            vertices = [int(a)-1 for a in f.readline().split()]
            face_vtx.append(vertices)
            connec[idx+1:idx+1+nside[i]] = vertices
            idx += nside[i] + 1

          # quick and dirty
          cell_type = [CellType.POLYGON for i in range(n_elt)]
          #grid = pv.UnstructuredGrid(connec, cell_type, points)

          grid = {
            "points": points,
            "cells" : connec,
            "face_vtx": face_vtx,
            "cell_type": cell_type
          }
          break
      instants[i_step] = grid

  # read scalars
  for name in scalars:
    padding = scalars[name]["padding"]
    fmt = "{:0%d}" % padding
    for i_step in range(n_steps):
      var_filename = path + scalars[name]["prefix"] + "." + fmt.format(start + incr*i_step)
      with open(var_filename, "r") as f:
        while True:
          line = f.readline()
          if line == "": break

          match = re.search("coordinates", line)
          if match is not None:
            n_vtx = len(instants[i_step]["points"])
            values = np.empty(n_vtx, dtype=float)
            for i in range(n_vtx):
              values[i] = float(f.readline())
            instants[i_step][name] = values
  return instants





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
      # print(filename)
      if not same_view:
        i_col = i_file%n_col
        i_row = i_file//n_col
        p.subplot(i_row, i_col)

      # load mesh file
      try:
        # instants = read_ensight(filename)
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


def parse_visu_cell(cell):
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
  return {
    "files" : files,
    "fields": fields,
    "styles": styles
  }


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
    params = parse_visu_cell(cell)

    # for l in cell_lines:
    #   a = [k.rstrip().lstrip() for k in l.split(":")]

    #   if len(a[0]) > 0:
    #     files.append(a[0])

    #     if len(a) > 1:
    #       fields.append(a[1])
    #     else:
    #       fields.append("")

    #     if len(a) > 2:
    #       styles.append(a[2])
    #     else:
    #       styles.append("")

    visu_n_files(files     = params["files"],
                 fields    = params["fields"],
                 styles    = params["styles"],
                 same_view = args.same_view,
                 lighting  = not args.no_lighting)


  """
  Visualize dynamic data sets
  """
  @cell_magic
  @magic_arguments.magic_arguments()
  @magic_arguments.argument('--same_view', '-sv',
                            help='Show all files in same view',
                            action="store_true")
  def visualize_dynamic(self, line='', cell=None):
    args = magic_arguments.parse_argstring(self.visualize, line)

    # parse cell content
    params = parse_visu_cell(cell)

    """
    Interactive visu (TODO: move to separate function)
    """
    from ipywidgets import interactive, IntSlider, HBox
    import matplotlib.pyplot as plt
    from matplotlib.colors import Colormap
    import numpy as np
    import copy

    # read data sets
    datasets = dict()
    for name in params["files"]:
      datasets[name] = read_ensight(name)

    n_files = len(params["files"])
    same_view = args.same_view or (n_files < 2)

    # get number of time steps
    n_steps = 0
    for name in datasets:
      n_steps = max(n_steps, len(datasets[name]))

    # compute global bbox over all time steps
    xymin =  1e30 * np.ones(2)
    xymax = -1e30 * np.ones(2)
    for name in datasets:
      for i_step in datasets[name]:
        coords = datasets[name][i_step]["points"]
        xymin = np.minimum(xymin, np.amin(coords[:,:2], axis=0))
        xymax = np.maximum(xymax, np.amax(coords[:,:2], axis=0))
    xyeps = 1e-2 * np.amax(xymax - xymin)
    xymin -= xyeps
    xymax += xyeps
    xmin, ymin = xymin
    xmax, ymax = xymax

    # compute scalar extrema over all time steps
    for i_dataset in range(n_files):
      name = params["files" ][i_dataset]
      vmin =  1e30
      vmax = -1e30
      field = params["fields"][i_dataset]
      for i_step in datasets[name]:
        if field in datasets[name][i_step]:
          values = datasets[name][i_step][field]
          vmin = min(vmin, np.amin(values))
          _vmax = np.amax(values)
          if _vmax < 1e6:
            vmax = max(vmax, _vmax)
      datasets[name]["vmin"] = vmin
      datasets[name]["vmax"] = vmax

    # TODO: link extrema?

    def visu(step):
      figsize = (15, 10)
      dpi = 72

      if same_view:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
      else:
        n_col = 2
        n_row = n_files//n_col
        fig, axes = plt.subplots(n_row, n_col, figsize=figsize, dpi=dpi)

      for i_dataset in range(n_files):
        name  = params["files" ][i_dataset]
        field = params["fields"][i_dataset]
        style = params["styles"][i_dataset]

        coords = datasets[name][step]["points"]
        connec = datasets[name][step]["face_vtx"]
        if field in datasets[name][step]:
          values = datasets[name][step][field]
        else:
          values = None

        if not same_view:
          i_col = i_dataset%n_col
          i_row = i_dataset//n_col
          if n_row == 1:
            ax = axes[i_col]
          else:
            ax = axes[i_row, i_col]

        if values is None:
          color = "0.5"
          lw    = 1.4
        else:
          #cmap = Colormap("coolwarm")
          cmap = copy.copy(plt.cm.get_cmap('coolwarm'))
          cmap.set_over("w")#'0.5')
          if style == "points":
            """ax.scatter(coords[:,0], coords[:,1],
                       c=values,
                       cmap=cmap,
                       zorder=3)"""
            status_name = field[:-1] + "_status"
            if status_name in datasets[name][step]:
              status = datasets[name][step][status_name]
              mapped   = np.where(np.isclose(status, 0))
              unmapped = np.where(np.isclose(status, 1))
              ax.plot(coords[mapped,  0], coords[mapped,  1], "ko", mfc="k", ms=5, zorder=3)
              ax.plot(coords[unmapped,0], coords[unmapped,1], "ko", mfc="w", ms=5, zorder=3)

          tpc = ax.tripcolor(coords[:,0], coords[:,1], connec,
                             values,
                             shading="gouraud",
                             cmap=cmap,
                             vmin=datasets[name]["vmin"],
                             vmax=datasets[name]["vmax"],
                             zorder=1)

          cbar = fig.colorbar(tpc,
                              ax=ax,
                              orientation="horizontal")
          # cbar.ax.tick_params(labelsize=14)
          cbar.set_label(label=field, size=14)
          color = "k"
          lw    = 0.7

        ax.triplot(coords[:,0], coords[:,1], connec, "-", color=color, lw=lw, zorder=2)

        ax.set_aspect("equal")
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.axis("off")

      plt.show()

    step = IntSlider(min=0, max=n_steps-1, step=1, value=0, continuous_update=True)
    ui = HBox([step])

    interactive_plot = interactive(visu, step=step)
    out = interactive_plot.children[-1]
    out.layout.height = '500px'
    display(ui, out)



def load_ipython_extension(ipython):
  ipython.register_magics(VisuMagics)
