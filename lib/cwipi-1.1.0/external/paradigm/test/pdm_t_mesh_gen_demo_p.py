#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import Pypdm.Pypdm as PDM


import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-n",      "--n_vtx_seg", default=3)
parser.add_argument("-n_part", "--n_part",    default=1)
parser.add_argument("-v",      "--visu",      action="store_true")

args = parser.parse_args()

n      = int(args.n_vtx_seg)
n_part = int(args.n_part)
visu   = args.visu

comm = MPI.COMM_WORLD

i_rank = comm.rank


"""
UV sphere
"""
dmn_capsule = PDM.sphere_surf_gen_nodal(comm,
                                        nu       = 2*n,
                                        nv       = n,
                                        x_center = 0,
                                        y_center = 0,
                                        z_center = 0,
                                        radius   = 1)

if visu:
  dmn_capsule.dump_vtk(PDM._PDM_GEOMETRY_KIND_SURFACIC,
                       "visu_sphere_")

"""
Icosphere
"""
dmn_capsule = PDM.sphere_surf_icosphere_gen_nodal(comm,
                                                  n,
                                                  x_center = 0,
                                                  y_center = 0,
                                                  z_center = 0,
                                                  radius   = 1)

if visu:
  dmn_capsule.dump_vtk(PDM._PDM_GEOMETRY_KIND_SURFACIC,
                       "visu_icosphere_")


"""
Icoball
"""
dmn_capsule = PDM.sphere_vol_icosphere_gen_nodal(comm,
                                                 n,
                                                 x_center = 0,
                                                 y_center = 0,
                                                 z_center = 0,
                                                 radius   = 1)

if visu:
  dmn_capsule.dump_vtk(PDM._PDM_GEOMETRY_KIND_VOLUMIC,
                       "visu_icoball_")

"""
Hollow ball
"""
dmn_capsule = PDM.sphere_vol_hollow_gen_nodal(comm,
                                              n,
                                              n_layer         = 2*n,
                                              x_center        = 0,
                                              y_center        = 0,
                                              z_center        = 0,
                                              radius_int      = 1,
                                              radius_ext      = 2,
                                              geometric_ratio = 1.1)

if visu:
  dmn_capsule.dump_vtk(PDM._PDM_GEOMETRY_KIND_VOLUMIC,
                       "visu_hollow_ball_")


"""
Dcube nodal
"""
dico = {
  # "_PDM_MESH_NODAL_POINT"    : PDM._PDM_MESH_NODAL_POINT   ,
  # "_PDM_MESH_NODAL_BAR2"     : PDM._PDM_MESH_NODAL_BAR2    ,
  "_PDM_MESH_NODAL_TRIA3"    : PDM._PDM_MESH_NODAL_TRIA3   ,
  "_PDM_MESH_NODAL_QUAD4"    : PDM._PDM_MESH_NODAL_QUAD4   ,
  "_PDM_MESH_NODAL_POLY_2D"  : PDM._PDM_MESH_NODAL_POLY_2D ,
  "_PDM_MESH_NODAL_TETRA4"   : PDM._PDM_MESH_NODAL_TETRA4  ,
  "_PDM_MESH_NODAL_PYRAMID5" : PDM._PDM_MESH_NODAL_PYRAMID5,
  "_PDM_MESH_NODAL_PRISM6"   : PDM._PDM_MESH_NODAL_PRISM6  ,
  "_PDM_MESH_NODAL_HEXA8"    : PDM._PDM_MESH_NODAL_HEXA8   ,
  # "_PDM_MESH_NODAL_POLY_3D"  : PDM._PDM_MESH_NODAL_POLY_3D ,
  }


for t in dico:
  dcube = PDM.DCubeNodalGenerator(nx     = n,
                                  ny     = n,
                                  nz     = n,
                                  length = 1,
                                  zero_x = 0,
                                  zero_y = 0,
                                  zero_z = 0,
                                  t_elmt = dico[t],
                                  order  = 1,
                                  comm   = comm)

  dcube.set_random_factor(1.)

  dcube.compute()

  dmn_capsule = dcube.get_dmesh_nodal()

  if t == "_PDM_MESH_NODAL_TRIA3" \
  or t == "_PDM_MESH_NODAL_QUAD4" \
  or t == "_PDM_MESH_NODAL_POLY_2D":
    geom_kind = PDM._PDM_GEOMETRY_KIND_SURFACIC
  else:
    geom_kind = PDM._PDM_GEOMETRY_KIND_VOLUMIC

  if visu:
    dmn_capsule.dump_vtk(geom_kind,
                         f"visu_dcube_{t[16:].lower()}")


"""
Ngon
"""
mesh = PDM.generate_mesh_rectangle_ngon(comm        = comm,
                                        elt_type    = PDM._PDM_MESH_NODAL_POLY_2D,
                                        xmin        = 0.,
                                        ymin        = 0.,
                                        zmin        = 0.,
                                        lengthx     = 1.,
                                        lengthy     = 1.,
                                        n_x         = n,
                                        n_y         = n,
                                        n_part      = n_part,
                                        part_method = PDM._PDM_SPLIT_DUAL_WITH_HILBERT,
                                        random_factor=1.)

if visu:
  PDM.writer_wrapper(comm,
                     "mesh_gen_demo",
                     "rectangle_ngon",
                     mesh["pvtx_coord"],
                     mesh["pvtx_ln_to_gn"],
                     mesh["pface_edge_idx"],
                     mesh["pface_vtx"],
                     mesh["pface_ln_to_gn"])

"""
Polyvol gen
"""
# ...
