#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import sys
import Pypdm.Pypdm as PDM
from Pypdm.Pypdm import npy_pdm_gnum_dtype

comm = MPI.COMM_WORLD

i_rank = comm.rank
n_rank = comm.size

assert(n_rank <= 2)

# fout = open("ex_part_to_part_{}.log".format(i_rank), "w")

if i_rank == 0:
  g_num1 = [
    np.array([1, 5]).astype(PDM.npy_pdm_gnum_dtype),
    np.array([6, 3, 7]).astype(PDM.npy_pdm_gnum_dtype)
  ]
  g_num2 = [
    np.array([2, 5, 3]).astype(PDM.npy_pdm_gnum_dtype),
    np.array([9, 8, 12, 11]).astype(PDM.npy_pdm_gnum_dtype)
  ]
  part1_to_part2 = [
    np.array([4, 3]).astype(PDM.npy_pdm_gnum_dtype),
    np.array([11, 1, 7]).astype(PDM.npy_pdm_gnum_dtype)
  ]

else:
  g_num1 = [
    np.array([2]).astype(PDM.npy_pdm_gnum_dtype),
    np.array([4, 8]).astype(PDM.npy_pdm_gnum_dtype)
  ]
  g_num2 = [
    np.array([1, 6, 4]).astype(PDM.npy_pdm_gnum_dtype),
    np.array([10, 7, 13]).astype(PDM.npy_pdm_gnum_dtype)
  ]
  part1_to_part2 = [
    np.array([5]).astype(PDM.npy_pdm_gnum_dtype),
    np.array([10, 2]).astype(PDM.npy_pdm_gnum_dtype)
  ]

part1_to_part2_idx = [
  np.array([i for i in range(len(g)+1)]).astype(np.intc)
 for g in g_num1]
# fout.write("part1_to_part2_idx = {}\n".format(part1_to_part2_idx))

# for i, g in enumerate(g_num1):
#   fout.write("part1 {}: gnums {}\n".format(i, g))
# for i, g in enumerate(g_num2):
#   fout.write("part2 {}: gnums {}\n".format(i, g))

# fout.write("\n\n=== part_to_part ===\n")


ptp = PDM.PartToPart(comm,
                     g_num1,
                     g_num2,
                     part1_to_part2_idx,
                     part1_to_part2)

ref_lnum2 = ptp.get_referenced_lnum2()
# fout.write("ref_lnum2 = {}\n".format(ref_lnum2))

come_from = ptp.get_gnum1_come_from()
# for i, g in enumerate(g_num2):
#   fout.write("part {}\n".format(i))
#   idx = come_from[i]["come_from_idx"]
#   for j, k in enumerate(ref_lnum2[i]):
#     cf = come_from[i]["come_from"][idx[j]:idx[j+1]]
#     fout.write("  gnum2 {} referenced by {}\n".format(g[k-1], cf))

stride = 3
# part1_stride = [
#   np.array([stride for i in g]).astype(np.intc) for g in g_num1
# ]
part1_stride = stride


is_interlaced = True#False

if is_interlaced:
  part1_data = [
    np.array([[i+(j+1)*0.1 for j in range(stride)] for i in g]).flatten() for g in g_num1
  ]
else:
  part1_data = [
    np.array([[i+(j+1)*0.1 for i in g] for j in range(stride)]).flatten() for g in g_num1
  ]

# fout.write("part1_data = {}\n".format(part1_data))
for i, p1d in enumerate(part1_data):
  # fout.write("part1 {}\n".format(i))
  for j, g in enumerate(g_num1[i]):
    if is_interlaced:
      data = p1d[stride*j:stride*(j+1)]
    else:
      data = p1d[j::len(g_num1[i])]
    # fout.write("  gnum1 {}: send {} to {}\n".format(g,
    #                                                 data,
    #                                                 part1_to_part2[i][part1_to_part2_idx[i][j]:part1_to_part2_idx[i][j+1]]))

request = ptp.iexch(PDM._PDM_MPI_COMM_KIND_P2P,
                    PDM._PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                    part1_data,
                    part1_stride=part1_stride,
                    interlaced_str=is_interlaced)

part2_stride, part2_data = ptp.wait(request)
# fout.write("part2_stride = {}\n".format(part2_stride))
# fout.write("part2_data   = {}\n".format(part2_data))
# fout.write("\n")

for i, g in enumerate(g_num2):
  # fout.write("part2 {}\n".format(i))
  idx = come_from[i]["come_from_idx"]
  for j, k in enumerate(ref_lnum2[i]):
    cf = come_from[i]["come_from"][idx[j]:idx[j+1]]
    if is_interlaced:
      data = part2_data[i][stride*idx[j]:stride*idx[j+1]]
    else:
      data = part2_data[i][idx[j]::len(ref_lnum2[i])]
    # fout.write("  gnum2 {} recv {} from {}\n".format(g[k-1],
    #                                                  data,
    #                                                  cf))



# fout.write("\n\n=== gnum_location ===\n")

gnum_loc = PDM.GlobalNumberingLocation(len(g_num2),
                                       len(part1_to_part2),
                                       comm)

for i, g in enumerate(g_num2):
  gnum_loc.gnum_location_elements_set(i,
                                      len(g),
                                      g)

for i, g in enumerate(part1_to_part2):
  gnum_loc.gnum_location_requested_elements_set(i,
                                                len(g),
                                                g)
gnum_loc.gnum_location_compute()

# for i, p1p2 in enumerate(part1_to_part2):
#   fout.write("part2 {}\n".format(i))
#   location_idx, location = gnum_loc.gnum_location_get(i)
#   for j, g in enumerate(p1p2):
#     fout.write("  gnum2 {}: location(s) : ".format(g))
#     for k in range(location_idx[j],location_idx[j+1],3):
#       fout.write("{} ".format(location[k:k+3]))
#     fout.write("\n")


# fout.close()


print("[{}] -- End".format(i_rank))
