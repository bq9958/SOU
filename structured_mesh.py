# -*- coding: utf-8 -*-
"""
Created on 04/05/2023
@author: 曾导SJTU
"""

import numpy as np
import math
from precision_data import fp

class StructuredMesh:
    def __init__(self):
        pass
    
    def CreateMesh(self, dim, ncx, ncy, ncz=1):
        # Number of spatial dimensions
        self.dim = dim

        # Number of cells in each direction
        self.ncx = ncx
        self.ncy = ncy
        self.ncz = ncz

        if self.dim == 2:
            self.ncz = 1

        # Number of nodes in each direction
        self.nx = self.ncx + 1
        self.ny = self.ncy + 1
        self.nz = self.ncz + 1

        print('nx, ny, nz = ', self.nx, self.ny, self.nz)
        print('ncx, ncy, ncz = ', self.ncx, self.ncy, self.ncz)

    def CreateCoordinates(self, xmin, xmax, ymin, ymax, zmin=fp(0.0), zmax=fp(1.0)):
        
        # Domain coordinates
        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.zmin, self.zmax = zmin, zmax

        if self.dim == 2:
            self.zmin = fp(0.0)
            self.zmax = fp(1.0)

        # Mesh coordinates
        self.x  = np.zeros(self.nx, dtype=fp)
        self.y  = np.zeros(self.ny, dtype=fp)
        self.z  = np.zeros(self.nz, dtype=fp)

        # Cell-centered mesh coordinates
        self.xc = np.zeros(self.ncx, dtype=fp)
        self.yc = np.zeros(self.ncy, dtype=fp)
        self.zc = np.zeros(self.ncz, dtype=fp)

        # Mesh generation
        self.dx = (self.xmax-self.xmin)/float(self.nx-1)
        for i in range(self.nx):
            self.x[i] = self.xmin + float(i)*self.dx
        self.dy = (self.ymax-self.ymin)/float(self.ny-1)
        for i in range(self.ny):
            self.y[i] = self.ymin + float(i)*self.dy
        self.dz = (self.zmax-self.zmin)/float(self.nz-1)
        for i in range(self.nz):
            self.z[i] = self.zmin + float(i)*self.dz

        for i in range(self.ncx):
            self.xc[i] = 0.5*(self.x[i]+self.x[i+1])
        for i in range(self.ncy):
            self.yc[i] = 0.5*(self.y[i]+self.y[i+1])
        for i in range(self.ncz):
            self.zc[i] = 0.5*(self.z[i]+self.z[i+1])

        print('dx, dy, dz = ', self.dx, self.dy, self.dz)
        print('bbox xmin, xmax = ', xmin, xmax)
        print('bbox ymin, ymax = ', ymin, ymax)
        print('bbox zmin, zmax = ', zmin, ymax)
        # print('x ', self.x)
        # print('y ', self.y)
        # print('z ', self.z)
        print('bbox Lx, Ly, Lz = ', xmax-xmin, ymax-ymin, zmax-zmin)

    def CreateFieldMeshData(self):
        
        # temperature field
        self.t = np.zeros((self.ncx, self.ncy, self.ncz), dtype=fp)
        
        # Old temperature field
        self.t0 = np.zeros((self.ncx, self.ncy, self.ncz), dtype=fp)

        # Intermediate velocity fields
        # Initialized as 0
        self.initial_uf = fp(0.0)
        self.uf = np.zeros((self.nx,  self.ncy, self.ncz), dtype=fp)
        self.initial_vf = fp(1.0)
        self.vf = np.zeros((self.ncx, self.ny,  self.ncz), dtype=fp)
        self.initial_wf = fp(0.0)
        self.wf = np.zeros((self.ncx, self.ncy, self.nz ), dtype=fp)
        
    def SetInitialT(self, T):
        self.t = T * np.ones((self.ncx, self.ncy, self.ncz), dtype=fp)

    def SetInitialUF(self, uf):
        self.initial_uf = uf
        self.uf = uf * np.ones((self.nx, self.ncy, self.ncz), dtype=fp)

    def GetInitialUF(self):
        return self.initial_uf

    def SetInitialVF(self, vf):
        self.uf = vf * np.ones((self.ncx, self.ny, self.ncz), dtype=fp)

    def SetInitialWF(self, wf):
        self.uf = wf * np.ones((self.ncx, self.ncy, self.nz), dtype=fp)
        
    def CreateCoeffMeshData(self):
        
        # Set number of coefficients
        if self.dim == 2:
            self.ncoef   = 6 # aP, aE, aW, aN, aS, bsrc
        elif self.dim == 3:
            self.ncoef = 8 # aP, aE, aW, aN, aS, aT, aB, bsrc      
        
        # Coefficient storage positions in four dimensional array
        self.id_aP = 0
        self.id_aE = 1
        self.id_aW = 2
        self.id_aN = 3
        self.id_aS = 4
        if self.dim == 3:
            self.id_aT = 5
            self.id_aB = 6
        self.id_bsrc = self.ncoef - 1        
        
        # temperature coefficient
        self.ct = np.zeros((self.ncx, self.ncy, self.ncz, self.ncoef  ), dtype=fp)
        
        self.r = np.zeros((self.nx, self.ny, self.nz), dtype = fp)

    def CreateSimulationData(self):
        
        # Equations
        self.ieqn = 0
        self.eqn_conduction = 0
        self.eqn_flow = 1
        self.eqn_conduction_flow = 2
        
        # Number of nonlinear iterations
        self.nsteps = 1
        
        # Time step, do not set it to be 0.0
        self.dt = fp(1.0)
        
        # Simulation loop control variables
        self.stop_sim = False
        
        # Convection scheme
        # 0: Upwind; 1: CD; 2: Power-law; 3: SOU;
        self.conv_scheme = 0

    def Set_nsteps(self, nsteps):
        self.nsteps = nsteps  

    def Set_dt(self, dt):
        self.dt = dt

    def Set_conv_scheme(self, conv_scheme):
        self.conv_scheme = conv_scheme
    
    def Set_limiter_scheme(self, limiter_scheme):
        self.limiter_scheme = limiter_scheme
        
    def CreateSolvingMethodData(self):
        
        # Number of linear solver iterations for temperature
        self.niter_t = 10
        
        # Relaxation factor
        self.relax_t = fp(0.75)
            
        # Convergence tolerances of the linear solver iterations
        self.res_t = fp(1.e-2)
        
        # Total number of linear iterations
        self.total_linsol_iters = 0
        
        # Convergence tolerances of the nonlinear iterations
        self.temp_tol = fp(1.e-6)                

        # L2-norm of residual of nonlinear iterations
        self.l2_curr  = fp(0.0)
        self.l2_max   = fp(-1.e20)
        self.l2_max_t = fp(-1.e20)
        self.l2_t     = fp(0.0)

    def Set_temp_solver_param(self, niter_t, relax_t, res_t, temp_tol):
        self.niter_t  = niter_t
        self.relax_t  = relax_t
        self.res_t    = res_t
        self.temp_tol = temp_tol


class Fluid:
    def __init__(self, ncx, ncy, ncz):
        
        self.ncx = ncx
        self.ncy = ncy
        self.ncz = ncz
        
        # density and viscosity fields
        self.dens = np.ones((ncx, ncy, ncz), dtype=fp)
        self.mu   = np.ones((ncx, ncy, ncz), dtype=fp) 

        # condution coefficient and specific heat
        self.con   = fp(0.0)
        self.spht  = fp(0.0)

        # heat source
        self.heat_src = fp(0.0)

    def SetInitialDenMu(self, dens, mu):
        # density and viscosity fields
        self.dens = dens * np.ones((self.ncx, self.ncy, self.ncz), dtype=fp)
        self.mu   = mu   * np.ones((self.ncx, self.ncy, self.ncz), dtype=fp)

    def Set_con_spht(self, con, spht):
        self.con   = con
        self.spht  = spht 
        

class BoundaryCondition: # None, Wall, Inlet, Outlet
    def __init__(self):
        self.type = 0

class BoundaryConditionTemp: 
    def __init__(self):
        self.temp_type = 0
        self.t = fp(0.0)
        self.heat_flux = fp(0.0)


class FluidBoundary:
    def __init__(self, dim):

        # Face index
        self.fid_e = 0
        self.fid_w = 1
        self.fid_n = 2
        self.fid_s = 3
        self.fid_t = 4
        self.fid_b = 5

        # Physical boundary condition types
        self.bc_none   = 0
        self.bc_wall   = 1
        self.bc_inlet  = 2
        self.bc_outlet = 3

        self.bcid_none   = 0 # internal edges/faces
        self.bcid_xmin   = 1
        self.bcid_xmax   = 2
        self.bcid_ymin   = 3
        self.bcid_ymax   = 4
        self.bcid_zmin   = 5
        self.bcid_zmax   = 6
        self.num_bcs = 7
        
        self.bcs = []
        self.bcs_temp = []
        for _ in range(self.num_bcs):
            self.bcs.append(BoundaryCondition())
            self.bcs_temp.append(BoundaryConditionTemp())
        
        # Values of physical temperature boundary conditions
        self.temp_east   = fp(0.0)
        self.temp_west   = fp(0.0)
        self.temp_north  = fp(0.0)
        self.temp_south  = fp(0.0)
        self.temp_top    = fp(0.0)
        self.temp_bottom = fp(0.0)
        
        # temperature boundary condition types
        self.temp_bc_constant = 0
        self.temp_bc_heatflux = 1                

    # 对于某个cell的每一个faces，设置它的boundary
    def CreateBoundaryOfCellFaces(self, case):
        
        # define local variables related to the case object
        dim  = case.dim  
        ncx  = case.ncx  
        ncy  = case.ncy  
        ncz  = case.ncz  
        xmin = case.xmin 
        xmax = case.xmax 
        ymin = case.ymin 
        ymax = case.ymax 
        zmin = case.zmin 
        zmax = case.zmax
        x    = case.x
        y    = case.y
        z    = case.z         
        
        # Boundary condition IDs
        self.bcid = np.zeros((ncx, ncy, ncz, 2*dim), dtype=int)
        
        # eps
        eps = fp(1.e-12)
        
        for k in range(ncz):
            for j in range(ncy):
                for i in range(ncx):
                    
                    # west face
                    x0 = x[i]
                    if abs(x0-xmin) < eps:
                        self.bcid[i,j,k,self.fid_w] = self.bcid_xmin
                    
                    # east face
                    x0 = x[i+1]
                    if abs(x0-xmax) < eps:
                        self.bcid[i,j,k,self.fid_e] = self.bcid_xmax
                                            
                    # south face
                    y0 = y[j]
                    if abs(y0-ymin) < eps:
                        self.bcid[i,j,k,self.fid_s] = self.bcid_ymin
             
                    # north face
                    y0 = y[j+1]
                    if abs(y0-ymax) < eps:
                        self.bcid[i,j,k,self.fid_n] = self.bcid_ymax
                        
                    if dim == 3:
                        # bottom face
                        z0 = z[k]
                        if abs(z0-zmin) < eps:
                            self.bcid[i,j,k,self.fid_b] = self.bcid_zmin
                
                        # top face
                        z0 = z[k+1]
                        if abs(z0-zmax) < eps:
                            self.bcid[i,j,k,self.fid_t] = self.bcid_zmax                                        


    def CreateBoundaryData(self, dim, # 2D
                           input_bc_xmin,
                           input_bc_xmax,
                           input_bc_ymin,
                           input_bc_ymax):

        # internal edges / faces
        id = self.bcid_none
        self.bcs[id].type = self.bc_none

        # xmin boundary conditions
        id = self.bcid_xmin
        buff = input_bc_xmin
        if buff == 'inlet':
            self.bcs[id].type = self.bc_inlet
        elif buff == 'outlet':
            self.bcs[id].type = self.bc_outlet
        elif buff == 'wall':
            self.bcs[id].type = self.bc_wall
            
        # xmax boundary conditions
        id = self.bcid_xmax
        buff = input_bc_xmax
        if buff == 'inlet':
            self.bcs[id].type = self.bc_inlet
        elif buff == 'outlet':
            self.bcs[id].type = self.bc_outlet
        elif buff == 'wall':
            self.bcs[id].type = self.bc_wall

        # ymin boundary conditions
        id = self.bcid_ymin
        buff = input_bc_ymin
        if buff == 'inlet':
            self.bcs[id].type = self.bc_inlet
        elif buff == 'outlet':
            self.bcs[id].type = self.bc_outlet
        elif buff == 'wall':
            self.bcs[id].type = self.bc_wall
            
        # ymax boundary conditions
        id = self.bcid_ymax
        buff = input_bc_ymax
        if buff == 'inlet':
            self.bcs[id].type = self.bc_inlet
        elif buff == 'outlet':
            self.bcs[id].type = self.bc_outlet
        elif buff == 'wall':
            self.bcs[id].type = self.bc_wall

        # 3D
        

    def CreateBoundaryDataTemp(self, dim,
                           input_bc_xmin, input_temp_type_xmin, input_heat_flux_xmin,
                           input_bc_xmax, input_temp_type_xmax, input_heat_flux_xmax,
                           input_bc_ymin, input_temp_type_ymin, input_heat_flux_ymin,
                           input_bc_ymax, input_temp_type_ymax, input_heat_flux_ymax):

        # internal edges / faces
        id = self.bcid_none
        self.bcs_temp[id].t = fp(0.0)
        self.bcs_temp[id].heat_flux = fp(0.0)
        self.bcs_temp[id].temp_type = self.bc_none

        # xmin boundary conditions
        id = self.bcid_xmin
        buff = input_temp_type_xmin
        if buff == 'constant':
            self.bcs_temp[id].temp_type = self.temp_bc_constant
            self.bcs_temp[id].t = input_heat_flux_xmin
        elif buff == 'heat_flux':
            self.bcs_temp[id].temp_type = self.temp_bc_heatflux
            self.bcs_temp[id].heat_flux = input_heat_flux_xmin
            
        # xmax boundary conditions
        id = self.bcid_xmax
        buff = input_temp_type_xmax
        if buff == 'constant':
            self.bcs_temp[id].temp_type = self.temp_bc_constant
            self.bcs_temp[id].t = input_heat_flux_xmax
        elif buff == 'heat_flux':
            self.bcs_temp[id].temp_type = self.temp_bc_heatflux
            self.bcs_temp[id].heat_flux = input_heat_flux_xmax

        # ymin boundary conditions
        id = self.bcid_ymin
        buff = input_temp_type_ymin
        if buff == 'constant':
            self.bcs_temp[id].temp_type = self.temp_bc_constant
            self.bcs_temp[id].t = input_heat_flux_ymin
        elif buff == 'heat_flux':
            self.bcs_temp[id].temp_type = self.temp_bc_heatflux
            self.bcs_temp[id].heat_flux = input_heat_flux_ymin
            
        # ymax boundary conditions
        id = self.bcid_ymax
        buff = input_temp_type_ymax
        if buff == 'constant':
            self.bcs_temp[id].temp_type = self.temp_bc_constant
            self.bcs_temp[id].t = input_heat_flux_ymax
        elif buff == 'heat_flux':
            self.bcs_temp[id].temp_type = self.temp_bc_heatflux
            self.bcs_temp[id].heat_flux = input_heat_flux_ymax

        # 3D
        
        
class PostProcessing:
    def __init__(self):
        
        # Output frequencies
        self.res_freq = 1
        self.out_freq = 1
                
        # Output file names
        self.linsol_fname    = "lin.res"
        self.nonlinsol_fname = "nonlin.res"
        self.vtk_fname_temp  = "post_temp.vtk"

        # Output file handles
        self.linsol_fid    = 101
        self.nonlinsol_fid = 102
        self.vtk_fid       = 103

    def Set_res_out_freq(self, res_freq, out_freq):
        self.res_freq = res_freq
        self.out_freq = out_freq
        
    def WriteVTKCollocated_temp(self, case, fluid):

        # Build the local data for np array
        nx = case.nx
        ny = case.ny
        nz = case.nz

        ncx = case.ncx
        ncy = case.ncy
        ncz = case.ncz
        
        x  = case.x
        y  = case.y
        z  = case.z        

        t  = case.t

        # Open VTK output file
        with open(self.vtk_fname_temp, 'w') as vtk_fid:

            # Write header
            vtk_fid.write('# vtk DataFile Version 3.0\n')
            vtk_fid.write('flash 3d grid and solution\n')
            vtk_fid.write('ASCII\n')
            vtk_fid.write('DATASET RECTILINEAR_GRID\n')

            # Write mesh information
            vtk_fid.write(f"DIMENSIONS {nx} {ny} {nz}\n")
            vtk_fid.write(f"X_COORDINATES {nx} float\n")
            vtk_fid.write(' '.join(str(i) for i in x) + '\n')
            vtk_fid.write(f"Y_COORDINATES {ny} float\n")
            vtk_fid.write(' '.join(str(i) for i in y) + '\n')
            vtk_fid.write(f"Z_COORDINATES {nz} float\n")
            vtk_fid.write(' '.join(str(i) for i in z) + '\n')

            # Write cell data
            ncell = (nx-1)*(ny-1)*(nz-1)
            vtk_fid.write(f"CELL_DATA {ncell}\n")
            
            vtk_fid.write('{:s}'.format("FIELD FieldData 1\n"))

            # Write temperature data
            vtk_fid.write(f"t 1 {ncell} float\n")
            t_arr = np.ravel(t[:, :, :], order='F')
            vtk_fid.write(' '.join(str(i) for i in t_arr) + '\n')

    def WriteVTKCollocated_temp_Pe_L(self, case, fluid):

        # Build the local data for np array
        nx = case.nx
        ny = case.ny
        nz = case.nz

        ncx = case.ncx
        ncy = case.ncy
        ncz = case.ncz
        
        x  = case.x
        y  = case.y
        z  = case.z        

        t  = case.t
        
        con = fluid.con
        if abs(con) > 10000:
            Pe_L = fp(0.0)
        else:
            Pe_L = case.GetInitialUF() / con
        
        print("Check Pe_L ", Pe_L)
        
        # # Open temperature output files
        if ncy==1:
            with open(f'temp_x_{Pe_L}.dat', 'w') as file3:
                # Write temperature data along x-line
                j = 0
                k = 0
                for i in range(nx-1):
                    file3.write(f"{(x[i]+x[i+1])*0.5} {t[i,j,k]}\n")

            with open(f'analytical_temp_x_{Pe_L}.dat', 'w') as file3:
                # Write temperature data along x-line
                j = 0
                k = 0
                if abs(Pe_L) < fp(1.e-3):
                    for i in range(nx-1):
                        file3.write(f"{(x[i]+x[i+1])*0.5} {(x[i]+x[i+1])*0.5}\n")
                else:
                    for i in range(nx-1):
                        xtemp = (x[i]+x[i+1])*0.5
                        result = (math.exp(Pe_L*xtemp)-1) / (math.exp(Pe_L)-1)
                        file3.write(f"{xtemp} {result}\n")                               
                
        if ncx==1:
            with open(f'temp_y_{Pe_L}.dat', 'w') as file3:
                # Write temperature data along x-line
                i = 0
                k = 0
                for j in range(ny-1):
                    file3.write(f"{(y[j]+y[j+1])*0.5} {t[i,j,k]}\n")

            with open(f'analytical_temp_y_{Pe_L}.dat', 'w') as file3:
                # Write temperature data along x-line
                i = 0
                k = 0
                if abs(Pe_L) < fp(1.e-3):
                    for j in range(ny-1):
                        file3.write(f"{(y[j]+y[j+1])*0.5} {(y[j]+y[j+1])*0.5}\n")
                else:
                    for j in range(ny-1):
                        ytemp = (y[j]+y[j+1])*0.5
                        result = (math.exp(Pe_L*ytemp)-1) / (math.exp(Pe_L)-1)
                        file3.write(f"{ytemp} {result}\n")

    def WriteVTKCollocated_temp_Pe_L_center(self, case, fluid):

        # Build the local data for np array
        nx = case.nx      

        t  = case.t
        
        con = fluid.con
        if abs(con) > 10000:
            Pe_L = fp(0.0)
        else:
            Pe_L = case.GetInitialUF() / con
            
        print("Check Pe_L ", Pe_L)
        
        # 0: Upwind; 1: CD; 2: Power-law; 3: SOU (to be implemented);
        conv_scheme = case.conv_scheme
        
        print("Check conv_scheme ", conv_scheme)
        
        if conv_scheme == 0:
            out_file = 'center_temp_x_upwind.dat'
        elif conv_scheme == 1:
            out_file = 'center_temp_x_center.dat'
        elif conv_scheme == 3:
            out_file = 'center_temp_x_SOU.dat'
        
        # Open temperature output files
        with open(out_file, 'a') as file3:
            # Write temperature data at center point
            i = int(nx/2)-1
            j = 0
            k = 0
            if abs(Pe_L) < fp(1.e-3):
                file3.write(f"{Pe_L} {t[i,j,k]} {0.5}\n")
            else:
                xtemp = 0.5
                result = (math.exp(Pe_L*xtemp)-1) / (math.exp(Pe_L)-1)
                file3.write(f"{Pe_L} {t[i,j,k]} {result}\n")  
