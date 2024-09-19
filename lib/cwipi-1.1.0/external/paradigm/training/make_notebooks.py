import os
#import jupytext
import argparse

def run(directory, recursive, verbose):
  for w in os.walk(directory):
    for x in w[2]:
      file_name, file_extension = os.path.splitext(x)
      if file_extension == ".md":
        if verbose:
          print("jupytext on ", w[0]+"/"+x)
        err = os.system("jupytext --to ipynb " + w[0] + "/" + x)
        #print(err)

    if not recursive: break

if __name__ == '__main__':
  parser = argparse.ArgumentParser()

  parser.add_argument("-d", "--dir",       type=str, default="./")
  parser.add_argument("-r", "--recursive", action="store_true")
  parser.add_argument("-v", "--verbose",   action="store_true")

  args = parser.parse_args()

  run(os.path.abspath(args.dir), args.recursive, args.verbose)

