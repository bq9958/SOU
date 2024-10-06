# -*- coding: utf-8 -*-
"""
Created on 04/05/2023
@author: 曾导SJTU
"""

import numpy as np
import time

from precision_data import fp
from structured_mesh import *
from collocated_segregated import *

# clock
clock_begin = time.perf_counter() # define clock_begin

# main codes
dim = 2

ncx = 25 # ncx, use odd number
ncy = 1 # ncy
ncz = 1 # ncz

# BCs on cell-centers
xmin = fp(-0.5/(ncx-1))
xmax = fp(1.0+0.5/(ncx-1)) # xmax, make sure xc[ncx-1] - xc[0] = 1.0
ymin = fp(0.0)
ymax = fp(1.0)
zmin = fp(0.0)
zmax = fp(1.0)

case = StructuredMesh()
case.CreateMesh(dim, ncx, ncy, ncz)
case.CreateCoordinates(xmin, xmax, ymin, ymax, zmin, zmax)
case.CreateFieldMeshData()

uf = fp(1.0) # uf
case.SetInitialUF(uf)

tref = fp(0.5) # initial condition
case.SetInitialT(tref)
        
case.CreateCoeffMeshData()

case.CreateSimulationData()

# 0: Upwind; 1: CD; 2: Power-law; 3: SOU (to be implemented);
conv_scheme = 1 # conv_scheme
case.Set_conv_scheme(conv_scheme)

limiter_scheme = "vanleer"
case.Set_limiter_scheme(limiter_scheme)

nsteps = 100
case.Set_nsteps(nsteps)
nsteps = case.nsteps

case.CreateSolvingMethodData()

niter_t = 100
relax_t = fp(0.75)
res_t = fp(0.1)
temp_tol = fp(1.e-6)
case.Set_temp_solver_param(niter_t, relax_t, res_t, temp_tol)

fluid = Fluid(ncx, ncy, ncz)
dens = fp(1.0)
mu   = fp(1.0)
fluid.SetInitialDenMu(dens, mu)

con  = fp(1000000) # control Pe_L
spht = fp(0.0)
fluid.Set_con_spht(con, spht)

fluidboundary = FluidBoundary(dim)

# 对于某个cell的每一个faces，设置它的boundary
fluidboundary.CreateBoundaryOfCellFaces(case)

fluidboundary.CreateBoundaryData(dim,
                                'wall',  # xmin
                                'wall',  # xmax
                                'wall',  # ymin
                                'wall' ) # ymax

# BC of temp also need to be changed if we wanna switch the simulation of convection-diffusion 
# from x to y direction
if ncy==1:
    fluidboundary.CreateBoundaryDataTemp(dim,
                                    'wall', 'constant', fp(0.0), # xmin
                                    'wall', 'constant', fp(1.0), # xmax
                                    'wall', 'heat_flux', fp(0.0), # ymin
                                    'wall', 'heat_flux', fp(0.0)) # ymax
elif ncx==1:
    fluidboundary.CreateBoundaryDataTemp(dim,
                                    'wall', 'heat_flux', fp(0.0), # xmin
                                    'wall', 'heat_flux', fp(0.0), # xmax
                                    'wall', 'constant', fp(0.0), # ymin
                                    'wall', 'constant', fp(1.0)) # ymax    

post = PostProcessing()
post.WriteVTKCollocated_temp(case, fluid)

res_freq = 50
out_freq = 1000
post.Set_res_out_freq(res_freq, out_freq)

nonlinsol_fname = post.nonlinsol_fname
linsol_fname = post.linsol_fname

# Open mass residual file
nonlinsol_fid = open(nonlinsol_fname, 'w')
nonlinsol_fid.write("#it, walltime, l2_t/l2_max_t\n")
nonlinsol_fid.close()

linsol_fid = open(linsol_fname, 'w')
linsol_fid.write("#it_nl, it, tot_it, norm, init, max, rel\n")
linsol_fid.close()

for it in range(1, nsteps+1):
    if (it % res_freq == 0 or it == 1 or it == nsteps):
        print('')
        print('----------------------------')
        print('Begin iter = ', it)
        print('----------------------------')
    
    case.stop_sim = collocated_segregated(it, case, fluid, fluidboundary, post)

    if (it % out_freq == 0 or it == nsteps or case.stop_sim):
        post.WriteVTKCollocated_temp(case, fluid)

    if (it == 1 or it % res_freq == 0 or it == nsteps or case.stop_sim):
        print("it, walltime, l2_t/l2_max_t ", it, time.perf_counter() - clock_begin, case.l2_t/case.l2_max_t)

        with open(nonlinsol_fname, 'a') as nonlinsol_fid:
            nonlinsol_fid.write("{} {} {}\n".format(it, time.perf_counter() - clock_begin, case.l2_t/case.l2_max_t))

    if case.stop_sim:
        break    

# post.WriteVTKCollocated_temp_Pe_L(case, fluid)
post.WriteVTKCollocated_temp_Pe_L_center(case, fluid)

elapsed_time = time.perf_counter() - clock_begin
print("elapsed time: ", elapsed_time)
print(" ")