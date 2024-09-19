import sys
import glob
import os
import json
from subprocess import run, CalledProcessError
from IPython.core import magic_arguments
from IPython.core.magic import cell_magic, line_magic, Magics, magics_class


headers = {
  "c"       : '#include <stdlib.h>\n#include <stdio.h>\n#include <string.h>\n#include <assert.h>\n#include <math.h>\n\n#include "pdm.h"\n#include "pdm_config.h"\n#include "pdm_mpi.h"\n#include "pdm_printf.h"\n#include "pdm_writer_priv.h"\n',
  "fortran" : "",
  "python"  : ""
}


language_extension = {
  "c"      : "c",
  "fortran": "F90",
  "python" : "py"
}

language_compiler_o = {
  "c"      : "/opt/tools/openmpi/4.0.5-gnu831/bin/mpicc",
  "fortran": "/opt/tools/openmpi/4.0.5-gnu831/bin/mpif90",
  "python" : None
}

language_compiler_e = {
  "c"      : "/opt/tools/openmpi/4.0.5-gnu831/bin/mpicxx",
  "fortran": "/opt/tools/openmpi/4.0.5-gnu831/bin/mpif90",
  "python" : None
}

user = "bandrieu"
pdm_dir   = f"/stck/{user}/workspace/paradigm/paradigm"
build_dir = f"/stck/{user}/workspace/trainings/build/paradigm"

language_linker_o = {
  "c"      : f"-DDEBUG_CLASSE -I{pdm_dir} -I{build_dir} -I{pdm_dir}/test -I{build_dir}/src -I{pdm_dir}/src/.. -I{pdm_dir}/src -I{pdm_dir}/src/pario -I{pdm_dir}/src/ppart -I{pdm_dir}/src/io -I{pdm_dir}/src/mpi_wrapper -I{pdm_dir}/src/ext_wrapper -I{pdm_dir}/src/mesh -I{pdm_dir}/src/meshgen -I{pdm_dir}/src/struct -I{pdm_dir}/src/util -I{pdm_dir}/src/gpu -I{pdm_dir}/src/adapt -I/opt/tools/scotch/6.0.9-idx32-gnu831-ompi405/include -I/opt/tools/parmetis/4.0.3-gnu831-ompi405/include -I/opt/tools/metis/5.1.0-gnu831/include -I{build_dir}/src/io -I{build_dir}/src/mpi_wrapper/mpi -I{pdm_dir}/src/mpi_wrapper/mpi -I{pdm_dir}/src/mpi_wrapper/mpi/.. -I{pdm_dir}/src/mpi_wrapper/mpi/../.. -std=gnu99 -fPIC -funsigned-char -pedantic -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused -Wfloat-equal  -Wno-unknown-pragmas -O0 -g -o",
  "fortran": f"-DDEBUG_CLASSE -I{pdm_dir} -I{build_dir} -I{pdm_dir}/test -I{build_dir}/src -I{pdm_dir}/src/.. -I{pdm_dir}/src -I{pdm_dir}/src/pario -I{pdm_dir}/src/ppart -I{pdm_dir}/src/io -I{pdm_dir}/src/mpi_wrapper -I{pdm_dir}/src/ext_wrapper -I{pdm_dir}/src/mesh -I{pdm_dir}/src/meshgen -I{pdm_dir}/src/struct -I{pdm_dir}/src/util -I{pdm_dir}/src/gpu -I{pdm_dir}/src/adapt -I/opt/tools/scotch/6.0.9-idx32-gnu831-ompi405/include -I/opt/tools/parmetis/4.0.3-gnu831-ompi405/include -I/opt/tools/metis/5.1.0-gnu831/include -I{build_dir}/src/io -I{build_dir}/src/mpi_wrapper/mpi -I{pdm_dir}/src/mpi_wrapper/mpi -I{pdm_dir}/src/mpi_wrapper/mpi/.. -I{pdm_dir}/src/mpi_wrapper/mpi/../.. -fallow-argument-mismatch -fPIC -cpp -Wall -std=gnu -Wno-unused-dummy-argument -Wno-maybe-uninitialized -O0 -g -fcheck=bounds -fbacktrace -o",
  "python" : None
}

language_linker_e1 = {
  "c"      : "-std=c++11 -fPIC -funsigned-char -W -Wall -Wshadow -Wpointer-arith -Wwrite-strings -Wunused -Wno-long-long -Wfloat-equal  -Wno-unknown-pragmas -O0 -g",
  "fortran": "-fallow-argument-mismatch -fPIC -cpp -Wall -std=gnu -Wno-unused-dummy-argument -Wno-maybe-uninitialized -O0 -g -fcheck=bounds -fbacktrace",
  "python" : None
}

language_linker_e2 = {
  "c"      : f"-L/opt/tools/scotch/6.0.9-idx32-gnu831-ompi405/lib  -Wl,-rpath,/opt/tools/parmetis/4.0.3-gnu831-ompi405/lib:/opt/tools/metis/5.1.0-gnu831/lib:/opt/tools/scotch/6.0.9-idx32-gnu831-ompi405/lib:::::::::::::::::::::::::::::::::::::::::::::::::::::::: -lm {build_dir}/src/libpdmf.a {build_dir}/src/libpdm.a {build_dir}/src/mpi_wrapper/mpi/libpdm_mpi.a {build_dir}/src/libpdmf.a {build_dir}/src/libpdm.a {build_dir}/src/mpi_wrapper/mpi/libpdm_mpi.a /opt/tools/parmetis/4.0.3-gnu831-ompi405/lib/libparmetis.so /opt/tools/metis/5.1.0-gnu831/lib/libmetis.so -lptscotch -lptscotcherr -lscotch -lscotcherr {build_dir}/src/io/libpdm_io.a -lm -lgfortran -lquadmath ",
  "fortran": f"-L/opt/tools/scotch/6.0.9-idx32-gnu831-ompi405/lib  -Wl,-rpath,/opt/tools/parmetis/4.0.3-gnu831-ompi405/lib:/opt/tools/metis/5.1.0-gnu831/lib:/opt/tools/scotch/6.0.9-idx32-gnu831-ompi405/lib:::::::::::::::::::::::::::::::::::::::::::::::::::::::: -lm {build_dir}/src/libpdmf.a {build_dir}/src/libpdm.a {build_dir}/src/mpi_wrapper/mpi/libpdm_mpi.a {build_dir}/src/libpdmf.a {build_dir}/src/libpdm.a {build_dir}/src/mpi_wrapper/mpi/libpdm_mpi.a /opt/tools/parmetis/4.0.3-gnu831-ompi405/lib/libparmetis.so /opt/tools/metis/5.1.0-gnu831/lib/libmetis.so -lptscotch -lptscotcherr -lscotch -lscotcherr {build_dir}/src/io/libpdm_io.a -lm -lstdc++ ",
  "python" : None
}

# temporary files (separate code blocks, not merged files and executables)
# will be placed here:
tmp_dir = "./tmp"



def print_out_err(proc):
  """
  Print output and error of a process
  (adapted from IDRIS trainings material)
  """
  stdout = proc.stdout
  stderr = proc.stderr
  if stderr:
    sys.stderr.write(stderr)
    sys.stderr.flush()
  if stdout:
    sys.stdout.write(stdout)
    sys.stdout.flush()



@magics_class
class CodeMagics(Magics):
  """
  Code block with integer identifier (to merge with other code blocks)
  """
  @cell_magic
  @magic_arguments.magic_arguments()
  @magic_arguments.argument('--id', '-i',
                            help='Block identifier (order)',
                            type=int,
                            default=0)
  @magic_arguments.argument('--prefix', '-p',
                            help='File name prefix',
                            type=str,
                            default="tmp")

  def code_block(self, line="", cell=None):
    args = magic_arguments.parse_argstring(self.code_block, line)

    extension = "tmp"
    str_id = "%3.3d" % args.id
    self.filename = f"{tmp_dir}/{args.prefix}_{str_id}.{extension}"

    if not os.path.exists(tmp_dir):
      os.makedirs(tmp_dir)

    # open source file and write cell content
    print(f"code block written in {self.filename}")
    with open(self.filename, "w+b") as out:
      out.write(bytes(cell, "utf-8"))
      out.flush()


  """
  Merge code blocks with given prefix into a new file
  The merged code can also be compiled and run in an MPI environment
  """
  @line_magic
  @magic_arguments.magic_arguments()
  @magic_arguments.argument('--prefix', '-p',
                            help='File name prefix',
                            type=str,
                            default=None)
  @magic_arguments.argument('--language', '-l',
                            help='Language',
                            type=str,
                            choices=["c", "python", "fortran"],
                            default="python")
  @magic_arguments.argument('--n_rank', '-n',
                            help='Number of MPI ranks',
                            type=int,
                            default=0)
  @magic_arguments.argument('--clear', '-c',
                            help='Clear generated files',
                            action="store_true")
  @magic_arguments.argument('--verbose', '-v',
                            help='Print merged source code',
                            action="store_true")
  def merge_code_blocks(self, line):
    # get current environment
    self.env = os.environ.copy() # ?

    args = magic_arguments.parse_argstring(self.merge_code_blocks, line)

    extension = language_extension[args.language]
    source_name = f"{tmp_dir}/{args.prefix}.{extension}"

    # get ordered list of code blocks
    tmp_files = glob.glob(f"{tmp_dir}/{args.prefix}_*.tmp")

    if len(tmp_files) == 0:
      print("No tmp files were found, make sure you run all code_block cells before this one!")
      return

    # merge into a single file
    with open(source_name, "w") as outfile:
      # writer langage header
      outfile.write(headers[args.language])
      # merge tmp files
      for fname in sorted(tmp_files):
        with open(fname) as infile:
          outfile.write(infile.read())

    if args.verbose:
      with open(source_name, "r") as f:
        print(f"merged :\n{f.read()}")

    run_code = args.n_rank > 0

    # Compile & run
    if run_code:
      # Compile
      if args.language == "python":
        exec_name = source_name
      else:
        exec_name = args.prefix
        o_name    = os.path.basename("{:s}.{:s}".format(args.prefix, "o"))

        command = []

        command.extend([language_compiler_o[args.language]])
        command.extend([language_linker_o[args.language]])
        command.extend([o_name, '-c', source_name])

        try:
          proc = run(" ".join(command),
                     capture_output=True,
                     shell=True,
                     check=True,
                     env=self.env,
                     encoding='utf8')
          print_out_err(proc)
        except CalledProcessError as e:
          sys.stderr.write(" ".join(e.cmd))
          sys.stderr.write(e.stderr)
          sys.stderr.write(e.stdout)
          return

        command = []

        command.extend([language_compiler_e[args.language]])
        command.extend([language_linker_e1[args.language]])
        command.extend([o_name, '-o', os.path.basename(exec_name)])
        command.extend([language_linker_e2[args.language]])

        try:
          proc = run(" ".join(command),
                     capture_output=True,
                     shell=True, # access to /d/whoami/workspace if True else only /stck/workspace
                     check=True,
                     env=self.env,
                     encoding='utf8')
          print_out_err(proc)
        except CalledProcessError as e:
          sys.stderr.write(" ".join(e.cmd))
          sys.stderr.write(e.stderr)
          sys.stderr.write(e.stdout)
          return

      # Run
      command = []
      command.extend(["mpirun"])
      command.extend(["-np"])
      command.extend(["%d" % args.n_rank])
      if args.language == "python":
        command.extend(["python3", "-u"])
      command.extend([exec_name])

      # sys.stdout.write(" ".join(command)+"\n")

      if os.path.isfile(exec_name):
        try:
          proc = run(command,
                     capture_output=True,
                     shell=False,
                     check=True,
                     env=self.env,
                     encoding='utf8')
          print_out_err(proc)
        except CalledProcessError as e:
          sys.stderr.write(" ".join(e.cmd))
          sys.stderr.write(e.stderr)
          sys.stderr.write(e.stdout)
          return
      else:
        sys.stderr.write(f"No executable was produced. It should be {exec_name}")

    # Clear generated files
    if args.clear:
      rm_files = tmp_files + [source_name]
      if args.language != "python":
        rm_files.append(o_name)
      if exec_name != source_name:
        rm_files += exec_name
      for fname in rm_files:
        if os.path.exists(fname):
          os.remove(fname)

def load_ipython_extension(ipython):
  ipython.register_magics(CodeMagics)
