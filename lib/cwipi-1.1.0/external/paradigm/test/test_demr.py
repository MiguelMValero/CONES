#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import sys
import Pypdm.Pypdm as PDM
from Pypdm.Pypdm import npy_pdm_gnum_dtype
import argparse

import sys
#sys.path.append("/stck/bandrieu/Public/adaptation/cavity_operator3_v3/")
# from mod_vtk import *


parser = argparse.ArgumentParser()

parser.add_argument("-n", "--n_vtx_seg", default=3)
parser.add_argument("-v", "--visu", action="store_true")

args = parser.parse_args()

n = int(args.n_vtx_seg)
visu = args.visu


comm = MPI.COMM_WORLD

i_rank = comm.rank
n_rank = comm.size


############################################################################
# On commence par générer un maillage sphérique "nodal" distribué
dmn_capsule = PDM.sphere_surf_icosphere_gen_nodal(comm,
                                                  n,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  1.)

# On partitionne ce maillage
n_part = 1
mpart = PDM.MultiPart(1,
                      np.array([n_part]).astype(np.intc),
                      0,
                      PDM._PDM_SPLIT_DUAL_WITH_HILBERT,
                      1,
                      np.ones(1).astype(np.double),
                      comm)

mpart.dmesh_nodal_set(0, dmn_capsule)

mpart.compute()
############################################################################


############################################################################
# On sélectionne un ensemble d'arêtes à partir desquelles on va générer des
# segments orthogonaux
extract_fraction = 0.25 # on sélectionne 25% des arêtes

# on récupère la connectivité arête->sommets des partitions
_, edge_vtx = mpart.connectivity_get(0, 0, PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX)

face_edge_idx, face_edge = mpart.connectivity_get(0, 0, PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE)

edge_face_idx, edge_face = mpart.connectivity_get(0, 0, PDM._PDM_CONNECTIVITY_TYPE_EDGE_FACE)

# part_infos = mpart.val_get(0, 0)

ghost_info_vtx  = mpart.ghost_information_get(0, 0)

# on récupère les numéros absolus des arêtes
edge_ln_to_gn = mpart.ln_to_gn_get(0, 0, PDM._PDM_MESH_ENTITY_EDGE)
n_edge = len(edge_ln_to_gn)

# on récupère les numéros absolus des sommets
vtx_ln_to_gn = mpart.ln_to_gn_get(0, 0, PDM._PDM_MESH_ENTITY_VTX)
# vtx_ln_to_gn = part_infos["np_vtx_ln_to_gn"]

# on récupère les coordonnées des sommets
vtx_coord = mpart.vtx_coord_get(0, 0)
n_vtx = len(vtx_ln_to_gn)

# vtx_coord = part_infos["np_vtx_coord"]

# on sélectionne les arêtes dont le numéro absolu est un multiple de 'step_edge'
step_edge = int(1./extract_fraction)
select_edge = []
for i in range(n_edge):
  if edge_ln_to_gn[i] % step_edge == 0:
    select_edge.append(i)


# Les arêtes sélectionnées sont extraites et redistribuées sur l'ensemble des processus
# pour générer les segments (fonctionnalité extract_part)

dim = 1 # les entités à extraire sont de dimension 1
extrp = PDM.ExtractPart(dim,
                        n_part,
                        n_part,
                        PDM._PDM_EXTRACT_PART_KIND_REEQUILIBRATE,#FROM_TARGET,
                        PDM._PDM_SPLIT_DUAL_WITH_HILBERT,
                        1,
                        comm)

# target_gnum = edge_ln_to_gn[select_edge]
# # pas nécessaire...
# # target_location = np.zeros(3*len(select_edge)).astype(np.intc)
# # for i, j in enumerate(select_edge):
# #   target_location[3*i  ] = i_rank
# #   target_location[3*i+2] = j

# target_location = None

# extrp.target_set(0, target_gnum, target_location)

extrp.selected_lnum_set(0, np.array(select_edge).astype(np.intc))

extrp.part_set(0,
               0, # n_cell
               0, # n_face
               n_edge,
               n_vtx,
               None, # cell_face_idx
               None, # cell_face
               None, # face_edge_idx
               None, # face_edge
               edge_vtx,
               None, # face_vtx_idx
               None, # face_vtx
               None, # cell_ln_to_gn
               None, # face_ln_to_gn
               edge_ln_to_gn,
               vtx_ln_to_gn,
               vtx_coord)

extrp.compute()

# on récupère les arêtes extraites et redistribuées
n_segment = extrp.n_entity_get(0, PDM._PDM_MESH_ENTITY_EDGE)

extrp_edge_vtx_idx, extrp_edge_vtx = extrp.connectivity_get(0, PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX)

# problème ici...
# segment_ln_to_gn = extrp.ln_to_gn_get(0, PDM._PDM_MESH_ENTITY_EDGE)
# print("segment_ln_to_gn = {}".format(segment_ln_to_gn))
# segment_parent_ln_to_gn = extrp.parent_ln_to_gn_get(0, PDM._PDM_MESH_ENTITY_EDGE)
# print("segment_parent_ln_to_gn = {}".format(segment_parent_ln_to_gn))

extrp_vtx_coord = extrp.vtx_coord_get(0)

segment_base_coord = np.zeros(3*n_segment)
segment_coord = np.zeros((2*n_segment, 3))
segment_connec = []
for i in range(n_segment):
  for j in range(2):
    for k in range(3):
      vtx_id = extrp_edge_vtx[2*i+j] - 1
      segment_base_coord[3*i+k] += 0.5*extrp_vtx_coord[3*vtx_id+k]
  segment_coord[2*i  ] = segment_base_coord[3*i:3*(i+1)]
  segment_coord[2*i+1] = 1.1*segment_coord[2*i]
  segment_connec.append([2*i, 2*i+1])

# on (re-)construit une numérotation absolue pour les segments
gen_gnum = PDM.GlobalNumbering(3, 1, 0, 1., comm)
gen_gnum.set_from_coords(0, segment_base_coord, np.ones(n_segment))
gen_gnum.compute()

segment_ln_to_gn = gen_gnum.get(0)
############################################################################



# A partir d'ici on a un maillage (notamment des arêtes) et des segments
# (en vision "partitionnée") répartis sur tous les processus


############################################################################
# On cherche maintenant à associer les segments avec les bonnes arêtes
# Pour cela, on construit 2 nuages de points:
# 1) la base des segments (cible)
# 2) le milieu des arêtes (source)
# Ensuite, pour chaque point du nuage 1 on cherche le point du nuage 2 le plus proche
n_closest_point = 1
clsp = PDM.ClosestPoints(comm, n_closest_point)

n_part_tgt = n_part
n_part_src = n_part
clsp.n_part_cloud_set(n_part_src, n_part_tgt)


# nuage source
edge_mid_coord = np.zeros(3*n_edge)
for i in range(n_edge):
  for j in range(2):
    vtx_id = edge_vtx[2*i+j] - 1
    edge_mid_coord[3*i:3*(i+1)] = edge_mid_coord[3*i:3*(i+1)] + 0.5*vtx_coord[3*vtx_id:3*(vtx_id+1)]

clsp.src_cloud_set(0, n_edge,    edge_mid_coord,     edge_ln_to_gn)

# nuage cible
clsp.tgt_cloud_set(0, n_segment, segment_base_coord, segment_ln_to_gn)

clsp.compute()

closest_src = clsp.points_get(0)
for i in range(n_segment):
  print("segment {} : edge {}, dist = {}".format(segment_ln_to_gn[i],
                                                 closest_src["closest_src_gnum"][i],
                                                 np.sqrt(closest_src["closest_src_distance"][i])))
############################################################################


# on récupère l'objet part_to_part pour faire des échanges entre les arêtes et les segments associés
ptp = clsp.part_to_part_get()

# ref_lnum2 = ptp.get_referenced_lnum2()
# print("ref_lnum2 = {}".format(ref_lnum2[0]))
# come_from = ptp.get_gnum1_come_from()
# print("come_from = {}".format(come_from[0]["come_from"]))

# on définit un champ basé sur les arêtes
edge_field = np.cos(10*((edge_mid_coord[0::3]+0.1)*(edge_mid_coord[1::3]+0.2)*(edge_mid_coord[2::3]+0.3)))

# on l'échange vers les segments
request = ptp.iexch(PDM._PDM_MPI_COMM_KIND_P2P,
                    PDM._PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                    [edge_field],
                    part1_stride=1,
                    interlaced_str=True)

part2_stride, part2_data = ptp.wait(request)

# print("edge_field = {}".format(edge_field))
# print("part2_data = {}".format(part2_data))


############################################################################
# Export for visu (ENSIGHT/VTK?)
# if visu:
#   faces = mpart.connectivity_get(0, 0, PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE)
#   face_edge_idx = faces["np_entity1_entity2_idx"]
#   face_edge     = faces["np_entity1_entity2"]

#   n_face = len(face_edge_idx)-1

#   face_vtx = PDM.compute_face_vtx_from_face_and_edge(face_edge_idx,
#                                                      face_edge,
#                                                      edge_vtx)


#   _edge_vtx = np.reshape(edge_vtx-1, (n_edge, 2), order="C")
#   _vtx_coord = np.reshape(vtx_coord, (n_vtx, 3), order="C")
#   vtk_write_std_elements("sphere_edges_rank%d.vtk" % i_rank,
#                          _vtx_coord,
#                          ELT_TYPE_EDGE,
#                          _edge_vtx,
#                          {"gnum" : edge_ln_to_gn,
#                          "field" : edge_field})

#   vtk_write_std_elements("sphere_segments_rank%d.vtk" % i_rank,
#                          segment_coord,
#                          ELT_TYPE_EDGE,
#                          segment_connec,
#                          {"gnum" : segment_ln_to_gn,
#                          "field" : part2_data[0]})


#   connec = []
#   for j in range(n_face):
#     start = face_edge_idx[j]
#     end   = face_edge_idx[j+1]
#     connec.append([k for k in face_vtx[start:end] - 1])

#   vtk_write_polydata("sphere_mesh_rank%d.vtk" % i_rank,
#                      _vtx_coord,
#                      connec)
############################################################################

print("[{}] End".format(i_rank))

