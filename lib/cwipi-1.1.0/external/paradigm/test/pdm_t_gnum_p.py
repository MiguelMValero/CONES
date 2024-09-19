#!/usr/bin/env python

def run_test():
  import numpy as np
  import argparse
  import mpi4py.MPI as MPI
  import Pypdm.Pypdm as PDM

  """
  Parse command line arguments
  """
  parser = argparse.ArgumentParser()

  parser.add_argument("-v",
                      "--verbose",
                      help="Verbose mode",
                      action="store_true")
  parser.add_argument("-m",
                      "--merge",
                      help="Merge quasi-coincident points",
                      action="store_true")
  parser.add_argument("-n_part",
                      "--n_part",
                      help="Number of partitions",
                      type=int,
                      default=1)
  parser.add_argument("-n",
                      "--n_pts",
                      help="Global number of points per partition",
                      type=int,
                      default=10)

  args = parser.parse_args()

  verbose = args.verbose
  merge   = 1 if args.merge else 0
  n_part  = int(args.n_part)
  n_pts   = int(args.n_pts)

  gnum = [None for i in range(n_part)]
  char_length = np.ones(n_pts) * 1e-3

  """
  Generate a random point cloud
  """
  n_rank = MPI.COMM_WORLD.size
  i_rank = MPI.COMM_WORLD.rank
  coords = []
  for i_part in range(n_part):
    cloud = PDM.part_point_cloud_gen_random(MPI.COMM_WORLD,
                                            n_rank*i_part + i_rank, # seed
                                            0,
                                            n_pts,
                                            0., 0., 0.,
                                            1., 1., 1.)

    coords.append(cloud["np_pts_coord"])

  """
  Generate a global numbering for the points from their coordinates
  """

  # First, create a GlobalNumbering instance and set some parameters
  gen_gnum = PDM.GlobalNumbering(3,    # dim
                                 n_part,
                                 merge,
                                 1e-3, # tolerance
                                 MPI.COMM_WORLD)

  # Then, provide the coordinates array for each partition
  # (`coords` is a list of numpy arrays of type double)
  # (`char_length` can be set to None if `merge` is disabled)
  for i_part in range(n_part):
    gen_gnum.set_from_coords(i_part,
                             coords[i_part],
                             char_length)

  # Once all partitions have been set, build the global numbering
  gen_gnum.compute()

  # Finally, retrieve the computed global id arrays
  for i_part in range(n_part):
    gnum[i_part] = gen_gnum.get(i_part)

    if verbose:
      print(f"part {i_part}: {gnum[i_part]}")


if __name__ == '__main__':
  run_test()
