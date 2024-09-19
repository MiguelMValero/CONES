import os, stat, sys
import itertools
import glob
import fnmatch
from pyparsing import pyparsing_common as ppc
from pyparsing import *


# def hooks_test(path):
#   """
#   """
#   # files = glob.glob(path + 'pdm_t_*.c', recursive=True)
#   # print(os.getcwd())
#   # print("path = ", path)
#   # print("files = ", files)
#   # # Filter only files
#   # files = [f for f in files if os.path.isfile(f)]
def hooks_test(path):
  """
  """
  file_list = []
  for path, folders, files in os.walk(path):
    for file in files:
      print(file)
      if fnmatch.fnmatch(file, 'pdm_t_*.c'):
        print(os.path.join(path, file))
        file_list.append(os.path.join(path, file))
  return file_list

filename = "/stck/bmaugars/dev/dev-Tools/paradigma/paradigm/test/pdm_t_dcube_nodal_gen.c"

param_tag = Literal('@@@param')
args_tag = Literal('@@@args')
def fparam(key, value):
  return param_tag.setResultsName("tag") + Literal("[") + key.setResultsName("key") + Literal("]") + Literal(":") + value

def fargs(key, value):
  return args_tag.setResultsName("tag") + Literal("[") + key.setResultsName("key") + Literal("]") + Literal(":") + value

def parse_one_file(filename):
  """
  """
  with open(filename, "r") as f:
    body = f.read()
  eol  = Literal( ';' )
  comment = cppStyleComment()

  ident = Regex("[a-zA-Z][a-zA-Z0-9_]*")

  numbers = Group(ppc.number + ZeroOrMore(Literal(",").suppress() + ppc.number)).setResultsName("value", listAllMatches=True)
  library = Regex("[a-zA-Z0-9_\\-]*")
  librairies = Group(library + ZeroOrMore(Literal(",").suppress() + library)).setResultsName("value", listAllMatches=True)

  parser = (fparam(Literal('n_proc'), numbers)          \
          | fparam(Literal('c'), numbers)               \
          | fparam(Literal('t'), numbers)               \
          | fparam(Literal('n'), numbers)               \
          | fparam(Literal('l'), numbers)               \
          | fparam(Literal('s'), numbers)               \
          | fparam(Literal('n_part'), numbers)          \
          | fparam(Literal('ext_type'), numbers)        \
          | fargs(Literal('part_kind'), librairies)     \
          | fargs(Literal('multipart'), librairies)     \
          | fargs(Literal('tree_kind'), librairies))

  infos = {}
  for (token, start, end) in parser.scanString(body):
    # print(f"token = {token.dump()}")
    # print(f"token.value = {token.value}")
    if token.key is not None:
      key = (token.key,token.tag)
      infos[key] = list(token.value[0])
  # print(f"infos = {infos}")

  # Treat n_proc
  parameters = {}
  key = ("n_proc", "@@@param")
  if(key not in infos):
    return parameters
  parameters[key[0]] = []


  for value in infos[key]:
    parameters[key[0]].append("-np {value}".format(value=value))

  for (key, tag) in filter(lambda k : k[0] != "n_proc", infos):
    print(f"key = {key}, tag = {tag}")
    parameters[key] = []
    if 'param' in tag:
      for value in infos[(key, tag)]:
        parameters[key].append("-{key} {value}".format(key=key, value=value))
    if 'args' in tag:
      for value in infos[(key, tag)]:
        parameters[key].append("{args}".format(args=value))
  # print(f"parameters = {parameters}")

  return parameters

# parameters = parse_one_file(filename)
def generate_shell(base_path_shell, directory, testname, parameters):
  shell_name = f"{base_path_shell}/run_tests_{testname}.sh"
  shell_run = open(shell_name, 'w')
  os.chmod(shell_name, stat.S_IRWXU | stat.S_IXGRP | stat.S_IRGRP)

  for n_proc in parameters["n_proc"]:
    cmd0 = "execute_test 'mpirun {n_proc} {directory}/{testname} ".format(testname=testname, directory=directory, n_proc=n_proc)
    lkeys = filter(lambda k : k != "n_proc", parameters)
    lparams = [parameters[key] for key in filter(lambda k : k != "n_proc", parameters)]
    for params in itertools.product(*lparams):
      # print("params = ", params)

      name_param = ''
      for param in params:
        name_param += "_"
        split_param = param.split(" ")
        if(len(split_param) == 2):
          name_param += '{0}={1}'.format(split_param[0][1::], split_param[1])
        else:
          name_param += '{0}'.format(split_param[0][1::])

      cmd  = f"class_name=\"{testname}{name_param}\"\n"
      cmd += f"test_name=\"{testname}_n={n_proc[4::]}{name_param}\"\n"
      cmd += f"test_n_proc=\"{n_proc[4::]}\"\n"
      cmd += cmd0 + " " + " ".join(params)
      cmd += "'"
      # print(f"cmd = {cmd}")
      shell_run.write(cmd+"\n")

  shell_run.close()

head, tail = os.path.split(filename)
# print(f"head = {head}, tail = {tail}")
testname, ext = os.path.splitext(tail)
# print(f"testname = {testname}, ext = {ext}")

# Paradigm
# test_files_pdm  = hooks_test("./test/")
# test_files_pdma = hooks_test("./extensions/paradigma/test/")

def parse_files_and_create_parameter_dict(base_path_shell, directory):
  """
  """
  test_files = hooks_test(directory)
  for filename in test_files:
    head, tail = os.path.split(filename)
    testname, ext = os.path.splitext(tail)
    parameters = parse_one_file(filename)
    if(len(parameters) > 0):
      generate_shell(base_path_shell, directory, testname, parameters)

base_path_shell = "script/prepare_ci/"
parse_files_and_create_parameter_dict(base_path_shell, "./test")
parse_files_and_create_parameter_dict(base_path_shell, "./extensions/paradigma/test")

