# -*- coding: utf-8 -*-
"""
Created on 04/05/2023
@author: 曾导SJTU
"""

import numpy as np
from precision_data import fp
from kernels import *

# For temperature
def conduction_coefs(case, dim, ncx, ncy, ncz, ncoef, dt, spht, con, heat_src, x, y, z, dens, t, uf, vf, wf, ct):

    # print('---------------')
    # print('conduction coefs')
    # print('---------------')

    # define local variables related to the case object
    id_aP = case.id_aP
    id_aE = case.id_aE
    id_aW = case.id_aW
    id_aN = case.id_aN
    id_aS = case.id_aS
    if dim == 3:
        id_aT = case.id_aT
        id_aB = case.id_aB
    id_bsrc = case.id_bsrc
    
    conv_scheme = case.conv_scheme

    idt = fp(1.0) / dt

    # Get area & volume
    dx = case.dx
    dy = case.dy
    dz = case.dz
    area_x, area_y, area_z, vol = cal_area_vol(dx, dy, dz)

    idx = fp(1.0) / dx
    idy = fp(1.0) / dy
    idz = fp(1.0) / dz
    
    for k in range(ncz):
        for j in range(ncy):
            for i in range(ncx):
                
                # Initialize coeffs
                aE = fp(0.0); aW = fp(0.0);
                aN = fp(0.0); aS = fp(0.0);
                aT = fp(0.0); aB = fp(0.0);
                aP = fp(0.0)
      
                # Initialize source & source coefficient to zero
                sC    = heat_src  # constant heat source
                sP    = fp(0.0)
                bsrc  = fp(0.0)

                # Calculate the coeffs as the Patankar's book
                # East Face
                if(i==ncx-1):
                    rho   = dens[i,j,k]
                else:
                    rho   = fp(0.5)*(dens[i,j,k]+dens[i+1,j,k])

                mul = con; mur = con
                ul    = uf[i+1,j,k]
                ur    = uf[i+1,j,k]
                aE    = a_nb(conv_scheme, area_x, idx, ul, ur, mul, mur, rho, -fp(1.0))

                # West Face
                if(i==0):
                    rho   = dens[i,j,k]
                else:
                    rho   = fp(0.5)*(dens[i,j,k]+dens[i-1,j,k])

                mul = con; mur = con
                ul    = uf[i,j,k]
                ur    = uf[i,j,k]
                aW    = a_nb(conv_scheme, area_x, idx, ul, ur, mul, mur, rho, fp(1.0))

                # North Face
                if(j==ncy-1):
                    rho   = dens[i,j,k]
                else:
                    rho   = fp(0.5)*(dens[i,j,k]+dens[i,j+1,k])

                mul = con; mur = con
                ul  = vf[i,j+1,k]
                ur  = vf[i,j+1,k]
                aN  = a_nb(conv_scheme, area_y, idy, ul, ur, mul, mur, rho, -fp(1.0))

                # South Face
                if(j==0):
                    rho   = dens[i,j,k]
                else:
                    rho   = fp(0.5)*(dens[i,j,k]+dens[i,j-1,k])

                mul = con; mur = con
                ul  = vf[i,j,k]
                ur  = vf[i,j,k]
                aS  = a_nb(conv_scheme, area_y, idy, ul, ur, mul, mur, rho, fp(1.0))

                if dim==3:

                    # Top Face
                    if(k==ncz-1):
                        rho   = dens[i,j,k]
                    else:
                        rho   = fp(0.5)*(dens[i,j,k]+dens[i,j,k+1])

                    mul = con; mur = con
                    ul  = wf[i,j,k+1]
                    ur  = wf[i,j,k+1]
                    aT  = a_nb(conv_scheme, area_z, idz, ul, ur, mul, mur, rho, -fp(1.0))

                    # Bottom Face
                    if(k==0):
                        rho   = dens[i,j,k]
                    else:
                        rho   = fp(0.5)*(dens[i,j,k]+dens[i,j,k-1])

                    mul = con; mur = con
                    ul  = wf[i,j,k]
                    ur  = wf[i,j,k]
                    aB  = a_nb(conv_scheme, area_z, idz, ul, ur, mul, mur, rho, fp(1.0))

                rho   = dens[i,j,k]
                aP0   = rho*spht*vol*idt
                aP    = aP + aE + aW + aN + aS + aT + aB + aP0 - sP*vol

                # Store coefficients as our reference book
                ct[i,j,k,id_aP]   = aP
                ct[i,j,k,id_aE]   = -aE
                ct[i,j,k,id_aW]   = -aW
                ct[i,j,k,id_aN]   = -aN
                ct[i,j,k,id_aS]   = -aS
                if dim==3:
                    ct[i,j,k,id_aT]  = -aT
                    ct[i,j,k,id_aB]  = -aB

                bsrc  = sC*vol + aP0*t[i,j,k]
                ct[i,j,k,id_bsrc] = bsrc
          
    return

def conduction_coef_bcs(case, fluidboundary, dim, ncx, ncy, ncz, ncoef, dt, con, x, y, z, dens, t, ct):

    # print('---------------')
    # print('conduction coef bcs')
    # print('---------------')

    # define local variables related to the case object
    id_aP = case.id_aP
    id_aE = case.id_aE
    id_aW = case.id_aW
    id_aN = case.id_aN
    id_aS = case.id_aS
    if dim == 3:
        id_aT = case.id_aT
        id_aB = case.id_aB
    id_bsrc = case.id_bsrc
    
    xc = case.xc
    yc = case.yc
    if dim==3:
        zc = case.zc

    # define local variables related to the fluidboundary object
    bcid  = fluidboundary.bcid

    fid_e = fluidboundary.fid_e
    fid_w = fluidboundary.fid_w
    fid_n = fluidboundary.fid_n
    fid_s = fluidboundary.fid_s
    if dim==3:
        fid_t = fluidboundary.fid_t
        fid_b = fluidboundary.fid_b
    
    bcs = fluidboundary.bcs
    bcs_temp = fluidboundary.bcs_temp
    
    bc_none   = fluidboundary.bc_none  
    bc_wall   = fluidboundary.bc_wall  
    bc_inlet  = fluidboundary.bc_inlet 
    bc_outlet = fluidboundary.bc_outlet
    
    temp_bc_constant = fluidboundary.temp_bc_constant
    temp_bc_heatflux = fluidboundary.temp_bc_heatflux

    # Temperature BC Coefs
    for k in range(ncz):
        for j in range(ncy):
            for i in range(ncx):
      
                bcid_e = bcid[i,j,k,fid_e]
                bcid_w = bcid[i,j,k,fid_w]
                bcid_n = bcid[i,j,k,fid_n]
                bcid_s = bcid[i,j,k,fid_s]
                # 3D
                
                bc_e = bcs[bcid_e].type
                bc_w = bcs[bcid_w].type
                bc_n = bcs[bcid_n].type
                bc_s = bcs[bcid_s].type
                
                ttype_e = bcs_temp[bcid_e].temp_type
                ttype_w = bcs_temp[bcid_w].temp_type
                ttype_n = bcs_temp[bcid_n].temp_type
                ttype_s = bcs_temp[bcid_s].temp_type
                
                bc_te = bcs_temp[bcid_e].t
                bc_tw = bcs_temp[bcid_w].t
                bc_tn = bcs_temp[bcid_n].t
                bc_ts = bcs_temp[bcid_s].t
                
                bc_fluxe = bcs_temp[bcid_e].heat_flux
                bc_fluxw = bcs_temp[bcid_w].heat_flux
                bc_fluxn = bcs_temp[bcid_n].heat_flux
                bc_fluxs = bcs_temp[bcid_s].heat_flux

                # East face BC
                if bc_e == bc_inlet:
                    ct[i,j,k,id_aP]   -= ct[i,j,k,id_aE]
                    ct[i,j,k,id_bsrc] -= 2.0 * ct[i,j,k,id_aE] * bc_te
                    ct[i,j,k,id_aE]   = 0.0
                elif bc_e == bc_outlet:
                    ct[i,j,k,id_aP]   += ct[i,j,k,id_aE]
                    ct[i,j,k,id_aE]   = 0.0
                elif bc_e == bc_wall:
                    if ttype_e == temp_bc_constant:
                        ct[i,j,k,id_aP]   -= ct[i,j,k,id_aE]
                        ct[i,j,k,id_bsrc] -= 2.0 * ct[i,j,k,id_aE] * bc_te
                        ct[i,j,k,id_aE]   = 0.0
                    elif ttype_e == temp_bc_heatflux:
                        ct[i,j,k,id_aP]   += ct[i,j,k,id_aE]
                        ct[i,j,k,id_bsrc] += ct[i,j,k,id_aE] * bc_fluxe * 2.0 * (x[i+1] - xc[i]) / con
                        ct[i,j,k,id_aE]   = 0.0

                # West face BC
                if bc_w == bc_inlet:
                    ct[i, j, k, id_aP] -= ct[i, j, k, id_aW]
                    ct[i, j, k, id_bsrc] -= 2.0 * ct[i, j, k, id_aW] * bc_tw
                    ct[i, j, k, id_aW] = 0.0
                elif bc_w == bc_outlet:
                    ct[i, j, k, id_aP] += ct[i, j, k, id_aW]
                    ct[i, j, k, id_aW] = 0.0
                elif bc_w == bc_wall:
                    if ttype_w == temp_bc_constant:
                        ct[i, j, k, id_aP] -= ct[i, j, k, id_aW]
                        ct[i, j, k, id_bsrc] -= 2.0 * ct[i, j, k, id_aW] * bc_tw
                        ct[i, j, k, id_aW] = 0.0
                    elif ttype_w == temp_bc_heatflux:
                        ct[i, j, k, id_aP] += ct[i, j, k, id_aW]
                        ct[i, j, k, id_bsrc] += ct[i, j, k, id_aW] * bc_fluxw * 2.0 * (x[i] - xc[i]) / con
                        ct[i, j, k, id_aW] = 0.0
                
                # North face BC
                if bc_n == bc_inlet:
                    ct[i, j, k, id_aP] -= ct[i, j, k, id_aN]
                    ct[i, j, k, id_bsrc] -= 2.0 * ct[i, j, k, id_aN] * bc_tn
                    ct[i, j, k, id_aN] = 0.0
                elif bc_n == bc_outlet:
                    ct[i, j, k, id_aP] += ct[i, j, k, id_aN]
                    ct[i, j, k, id_aN] = 0.0
                elif bc_n == bc_wall:
                    if ttype_n == temp_bc_constant:
                        ct[i, j, k, id_aP] -= ct[i, j, k, id_aN]
                        ct[i, j, k, id_bsrc] -= 2.0 * ct[i, j, k, id_aN] * bc_tn
                        ct[i, j, k, id_aN] = 0.0
                    elif ttype_n == temp_bc_heatflux:
                        ct[i, j, k, id_aP] += ct[i, j, k, id_aN]
                        ct[i, j, k, id_bsrc] += ct[i, j, k, id_aN] * bc_fluxn * 2.0 * (y[j+1] - yc[j]) / con
                        ct[i, j, k, id_aN] = 0.0

                # South face BC
                if bc_s == bc_inlet:
                    ct[i,j,k,id_aP]   -= ct[i,j,k,id_aS]
                    ct[i,j,k,id_bsrc] -= 2.0*ct[i,j,k,id_aS]*bc_ts
                    ct[i,j,k,id_aS]    = 0.0
                elif bc_s == bc_outlet:
                    ct[i,j,k,id_aP]   += ct[i,j,k,id_aS]
                    ct[i,j,k,id_aS]    = 0.0
                elif bc_s == bc_wall:
                    if ttype_s == temp_bc_constant:
                        ct[i,j,k,id_aP]   -= ct[i,j,k,id_aS]
                        ct[i,j,k,id_bsrc] -= 2.0*ct[i,j,k,id_aS]*bc_ts
                        ct[i,j,k,id_aS]    = 0.0
                    elif ttype_s == temp_bc_heatflux:
                        ct[i,j,k,id_aP]   += ct[i,j,k,id_aS]
                        ct[i,j,k,id_bsrc] += ct[i,j,k,id_aS]*bc_fluxs*2.0*(y[j]-yc[j])/con
                        ct[i,j,k,id_aS]    = 0.0

                # 3D

    # Forced cell-centered boundary conditions for 红宝书的例子，例子见书上Fig. 11.7
    # for k in range(ncz):
    #     for j in range(ncy):
    #         for i in range(ncx):
    #             if i==0:
    #                 ct[i,j,k,id_aP]   = fp(1.0)
    #                 ct[i,j,k,id_aE]   = fp(0.0)
    #                 ct[i,j,k,id_aW]   = fp(0.0)
    #                 ct[i,j,k,id_aN]   = fp(0.0)
    #                 ct[i,j,k,id_aS]   = fp(0.0)
    #                 ct[i,j,k,id_bsrc] = bc_tw
    #             elif i==ncx-1:
    #                 ct[i,j,k,id_aP]   = fp(1.0)
    #                 ct[i,j,k,id_aE]   = fp(0.0)
    #                 ct[i,j,k,id_aW]   = fp(0.0)
    #                 ct[i,j,k,id_aN]   = fp(0.0)
    #                 ct[i,j,k,id_aS]   = fp(0.0)
    #                 ct[i,j,k,id_bsrc] = bc_te

    return

def SOU_src(case, fluidboundary, dim, ncx, ncy, ncz, t, uf, vf, wf, dens, ct):
    
    limiter_scheme = case.limiter_scheme
    
    id_bsrc = case.id_bsrc
    
    # Get area & volume
    dx = case.dx
    dy = case.dy
    dz = case.dz
    area_x, area_y, area_z, vol = cal_area_vol(dx, dy, dz)
    
    xc = case.xc
    yc = case.yc
    if dim==3:
        zc = case.zc

    # define local variables related to the fluidboundary object
    bcid  = fluidboundary.bcid

    fid_e = fluidboundary.fid_e
    fid_w = fluidboundary.fid_w
    fid_n = fluidboundary.fid_n
    fid_s = fluidboundary.fid_s
    if dim==3:
        fid_t = fluidboundary.fid_t
        fid_b = fluidboundary.fid_b
    
    bcs = fluidboundary.bcs
    bcs_temp = fluidboundary.bcs_temp
    
    bc_none   = fluidboundary.bc_none  
    bc_wall   = fluidboundary.bc_wall  
    bc_inlet  = fluidboundary.bc_inlet 
    bc_outlet = fluidboundary.bc_outlet
    
    temp_bc_constant = fluidboundary.temp_bc_constant
    temp_bc_heatflux = fluidboundary.temp_bc_heatflux
    
    for k in range(ncz):
            for j in range(ncy):
                for i in range(ncx):
                                        
                    src = fp(0.0)
                    
                    # phiuu | phiu | p | phid | phidd
                    # u > 0 toward east; v > 0 toward north
                    # therfore, west and south are upstreams; east and north are downstreams
                    t_uu = fp(0.)
                    t_u  = fp(0.)
                    t_p  = t[i,j,k]
                    t_d  = fp(0.)
                    t_dd = fp(0.)
                    
                    # face vel up/downstream
                    uf_u = uf[i-1,j,k]
                    uf_d = uf[i,j,k]
                    vf_u = vf[i,j-1,k]
                    vf_d = vf[i,j,k]
                    
                    current_cell_boundary = False
                    
                    rho = fp(0.)
                    
                    bcid_e = bcid[i,j,k,fid_e]
                    bcid_w = bcid[i,j,k,fid_w]
                    bcid_n = bcid[i,j,k,fid_n]
                    bcid_s = bcid[i,j,k,fid_s]
                    
                    bc_e = bcs[bcid_e].type
                    bc_w = bcs[bcid_w].type
                    bc_n = bcs[bcid_n].type
                    bc_s = bcs[bcid_s].type
                    
                    ttype_e = bcs_temp[bcid_e].temp_type
                    ttype_w = bcs_temp[bcid_w].temp_type
                    ttype_n = bcs_temp[bcid_n].temp_type
                    ttype_s = bcs_temp[bcid_s].temp_type
                    
                    # East&west faces
                                        
                    # if the east face of the current cell is boundary
                    if bc_e == bc_inlet:
                        #logging.debug("[East] Current cell is boundary")
                        if ttype_e == temp_bc_heatflux:
                            bc_fluxe = bcs_temp[bcid_e].heat_flux
                            t_d  = t_p + bc_fluxe*dx
                        elif ttype_e == temp_bc_constant:
                            bc_te = bcs_temp[bcid_e].t
                            t_d = fp(2.0)*bc_te-t_p 
                        t_dd = fp(2.0)*t_d-t_p
                        rho = dens[i,j,k]
                        current_cell_boundary = True
                    elif bc_e == bc_wall:
                        #logging.debug("[East] Current cell is boundary")
                        if ttype_e == temp_bc_heatflux:
                            bc_fluxe = bcs_temp[bcid_e].heat_flux
                            t_d  = t_p + bc_fluxe*dx
                        elif ttype_e == temp_bc_constant:
                            bc_te = bcs_temp[bcid_e].t
                            t_d = fp(2.0)*bc_te-t_p 
                        t_dd = fp(2.0)*t_d-t_p
                        rho = dens[i,j,k]
                        current_cell_boundary = True
                    elif bc_e == bc_outlet:
                        #logging.debug("[East] Current cell is boundary")
                        if ttype_e == temp_bc_heatflux:
                            bc_fluxe = bcs_temp[bcid_e].heat_flux
                            t_d  = t_p + bc_fluxe*dx
                        elif ttype_e == temp_bc_constant:
                            bc_te = bcs_temp[bcid_e].t
                            t_d = fp(2.0)*t_p-t_u 
                        t_dd = fp(2.0)*t_d-t_p
                        rho = dens[i,j,k]
                        current_cell_boundary = True
                    else:
                        t_d = t[i+1,k,j]
                        rho = fp(0.5) * (dens[i,j,k] + dens[i+1,j,k])
                    
                    # if the east face of the east cell is boundary
                    if current_cell_boundary == False:
                        #logging.debug("[East] Current cell is not boundary")
                        if bcs[bcid[i+1,j,k,fid_e]].type == bc_inlet:
                            #logging.debug("[East] Neighbor is boundary")
                            bcid_e = bcid[i+1,j,k,fid_e]
                            bc_te = bcs_temp[bcid_e].t
                            ttype_e = bcs_temp[bcid_e].temp_type
                            if ttype_e == temp_bc_heatflux:
                                bc_fluxe = bcs_temp[bcid_e].heat_flux
                                t_dd = t_d + bc_fluxe*dx
                            elif ttype_e == temp_bc_constant:
                                t_dd = fp(2.0)*bc_te-t_d
                            rho = fp(0.5) * (dens[i,j,k] + dens[i+1,j,k])
                        elif bcs[bcid[i+1,j,k,fid_e]].type == bc_wall:
                            #logging.debug("[East] Neighbor is boundary")
                            bcid_e = bcid[i+1,j,k,fid_e]
                            bc_te = bcs_temp[bcid_e].t
                            ttype_e = bcs_temp[bcid_e].temp_type
                            if ttype_e == temp_bc_heatflux:
                                bc_fluxe = bcs_temp[bcid_e].heat_flux
                                t_dd = t_d + bc_fluxe*dx
                            elif ttype_e == temp_bc_constant:
                                t_dd = fp(2.0)*bc_te-t_d
                            rho = fp(0.5) * (dens[i,j,k] + dens[i+1,j,k])
                        elif bcs[bcid[i+1,j,k,fid_e]].type == bc_outlet:
                            #logging.debug("[East] Neighbor is boundary")
                            bcid_e = bcid[i+1,j,k,fid_e]
                            bc_te = bcs_temp[bcid_e].t
                            ttype_e = bcs_temp[bcid_e].temp_type
                            if ttype_e == temp_bc_heatflux:
                                bc_fluxe = bcs_temp[bcid_e].heat_flux
                                t_dd = t_d + bc_fluxe*dx
                            elif ttype_e == temp_bc_constant:
                                t_dd = fp(2.0)*bc_te-t_d
                            rho = fp(0.5) * (dens[i,j,k] + dens[i+1,j,k])
                        else:
                            #logging.debug("[East] Neighbor is not boundary")
                            t_dd = t[i+2,k,j]
                            rho = fp(0.5) * (dens[i,j,k] + dens[i+1,j,k])
                        
                    f_d = rho*uf_d*area_x   # flux downstream
                    current_cell_boundary = False
                    
                    # if the west face of the current cell is boundary
                    if bc_w == bc_inlet:
                        #logging.debug("[West] Current cell is boundary")
                        if ttype_w == temp_bc_heatflux:
                            bc_fluxw = bcs_temp[bcid_w].heat_flux
                            t_u  = t_p + bc_fluxw*dx
                        elif ttype_w == temp_bc_constant:
                            bc_tw = bcs_temp[bcid_w].t
                            t_u = fp(2.0)*bc_tw-t_p 
                        t_uu = fp(2.0)*t_u-t_p
                        rho = dens[i,j,k]
                        current_cell_boundary = True
                    elif bc_w == bc_wall:
                        #logging.debug("[West] Current cell is boundary")
                        if ttype_w == temp_bc_heatflux:
                            bc_fluxw = bcs_temp[bcid_w].heat_flux
                            t_u  = t_p + bc_fluxw*dx
                        elif ttype_w == temp_bc_constant:
                            bc_tw = bcs_temp[bcid_w].t
                            t_u = fp(2.0)*bc_tw-t_p 
                        t_uu = fp(2.0)*t_u-t_p
                        rho = dens[i,j,k]
                        current_cell_boundary = True
                    elif bc_w == bc_outlet:
                        #logging.debug("[West] Current cell is boundary")
                        if ttype_w == temp_bc_heatflux:
                            bc_fluxw = bcs_temp[bcid_w].heat_flux
                            t_u  = t_p + bc_fluxw*dx
                        elif ttype_w == temp_bc_constant:
                            bc_tw = bcs_temp[bcid_w].t
                            t_u = fp(2.0)*bc_tw-t_p 
                        t_uu = fp(2.0)*t_u-t_p
                        rho = dens[i,j,k]
                        current_cell_boundary = True
                    else:
                        t_u = t[i-1,k,j]
                        rho = fp(0.5) * (dens[i-1,j,k] + dens[i,j,k])
                    
                    # if the west face of the west cell is boundary
                    if current_cell_boundary == False:
                        #logging.debug("[West] Current cell is not boundary")
                        if bcs[bcid[i-1,j,k,fid_e]].type == bc_inlet:
                            bcid_w = bcid[i-1,j,k,fid_e]
                            bc_tw = bcs_temp[bcid_w].t
                            ttype_w = bcs_temp[bcid_w].temp_type
                            if ttype_w == temp_bc_heatflux:
                                bc_fluxw = bcs_temp[bcid_w].heat_flux
                                t_uu = t_u + bc_fluxw*dx
                            elif ttype_w == temp_bc_constant:
                                t_uu = fp(2.0)*bc_tw-t_u
                            rho = fp(0.5) * (dens[i,j,k] + dens[i-1,j,k])
                        elif bcs[bcid[i-1,j,k,fid_w]].type == bc_wall:
                            #logging.debug("[West] Neighbor is boundary")
                            bcid_w = bcid[i-1,j,k,fid_e]
                            bc_tw = bcs_temp[bcid_w].t
                            ttype_w = bcs_temp[bcid_w].temp_type
                            if ttype_w == temp_bc_heatflux:
                                bc_fluxw = bcs_temp[bcid_w].heat_flux
                                t_uu = t_u + bc_fluxw*dx
                            elif ttype_w == temp_bc_constant:
                                t_uu = fp(2.0)*bc_tw-t_u
                            rho = fp(0.5) * (dens[i,j,k] + dens[i-1,j,k])
                        elif bcs[bcid[i-1,j,k,fid_w]].type == bc_outlet:
                            #logging.debug("[West] Neighbor is boundary")
                            bcid_w = bcid[i-1,j,k,fid_e]
                            bc_tw = bcs_temp[bcid_w].t
                            ttype_w = bcs_temp[bcid_w].temp_type
                            if ttype_w == temp_bc_heatflux:
                                bc_fluxw = bcs_temp[bcid_w].heat_flux
                                t_uu = t_u + bc_fluxw*dx
                            elif ttype_w == temp_bc_constant:
                                t_uu = fp(2.0)*bc_tw-t_u
                            rho = fp(0.5) * (dens[i,j,k] + dens[i-1,j,k])
                        else:
                            #logging.debug("[West] Neighbor is not boundary")
                            t_uu = t[i-2,j,k]
                            rho = fp(0.5) * (dens[i,j,k] + dens[i-1,j,k])
                        
                        
                    f_u = rho*uf_u*area_x   # flux upstream
                    current_cell_boundary = False
                    
                    au = fp(0.); ad = fp(0.)
                    if f_u > fp(0.):
                        au = fp(1.)
                    if f_d > fp(0.):
                        ad = fp(1.)
                    
                    rdm = (t_dd-t_d)/(t_d-t_p+1.e-12)    # rdm : r_downstream^minus
                    rdp = (t_p-t_u)/(t_d-t_p+1.e-12)     # rdp : r_downstream^plus
                    rum = (t_d-t_p)/(t_p-t_u+1.e-12)
                    rup = (t_u-t_uu)/(t_p-t_u+1.e-12)
                    
                    if(t_dd==t_d) : rdm = fp(0.)
                    if(t_p==t_u) : rdp = fp(0.)
                    if(t_d==t_p) : rum = fp(0.)
                    if(t_u==t_uu) : rup = fp(0.)
                    
                    src = (fp(0.5)*f_d*((fp(1.)-ad)*limiter(rdm,limiter_scheme)-ad*limiter(rdp,limiter_scheme))*(t_d-t_p) 
                                - fp(0.5)*f_u*((fp(1.)-au)*limiter(rum,limiter_scheme)-au*limiter(rup,limiter_scheme))*(t_p-t_u))
                    
                    # North&east faces
                    # face vel down/upstream
                        
                    # if the north face of the current cell is boundary
                    if bc_n == bc_inlet:
                        #logging.debug("[North] Current cell is boundary")
                        if ttype_n == temp_bc_heatflux:
                            bc_fluxn = bcs_temp[bcid_n].heat_flux
                            t_d  = t_p + bc_fluxn*dy
                        elif ttype_n == temp_bc_constant:
                            bc_tn = bcs_temp[bcid_n].t
                            t_d = fp(2.0)*bc_tn-t_p 
                        t_dd = fp(2.0)*t_d-t_p
                        rho = dens[i,j,k]
                        current_cell_boundary = True
                    elif bc_n == bc_wall:
                        #logging.debug("[North] Current cell is boundary")
                        if ttype_n == temp_bc_heatflux:
                            bc_fluxn = bcs_temp[bcid_n].heat_flux
                            t_d  = t_p + bc_fluxn*dy
                        elif ttype_n == temp_bc_constant:
                            bc_tn = bcs_temp[bcid_n].t
                            t_d = fp(2.0)*bc_tn-t_p 
                        t_dd = fp(2.0)*t_d-t_p
                        rho = dens[i,j,k]
                        current_cell_boundary = True
                    elif bc_n == bc_outlet:
                        #logging.debug("[North] Current cell is boundary")
                        if ttype_n == temp_bc_heatflux:
                            bc_fluxn = bcs_temp[bcid_n].heat_flux
                            t_d  = t_p + bc_fluxn*dy
                        elif ttype_n == temp_bc_constant:
                            bc_tn = bcs_temp[bcid_n].t
                            t_d = fp(2.0)*bc_tn-t_p 
                        t_dd = fp(2.0)*t_d-t_p
                        rho = dens[i,j,k]
                        current_cell_boundary = True
                    else:
                        t_d = t[i,k+1,j]
                        rho = fp(0.5) * (dens[i,j,k] + dens[i,j+1,k])
                    
                    # if the north face of the north cell is boundary
                    if current_cell_boundary == False:
                        if bcs[bcid[i,j+1,k,fid_n]].type == bc_inlet:
                            bcid_n = bcid[i,j+1,k,fid_n]
                            bc_tn = bcs_temp[bcid_n].t
                            ttype_n = bcs_temp[bcid_n].temp_type
                            if ttype_n == temp_bc_heatflux:
                                bc_fluxn = bcs_temp[bcid_n].heat_flux
                                t_dd = t_d + bc_fluxn*dy
                            elif ttype_n == temp_bc_constant:
                                t_dd = fp(2.0)*bc_tn-t_d
                            rho = fp(0.5) * (dens[i,j,k] + dens[i,j+1,k])
                        elif bcs[bcid[i,j+1,k,fid_n]].type == bc_wall:
                            bcid_n = bcid[i,j+1,k,fid_n]
                            bc_tn = bcs_temp[bcid_n].t
                            ttype_n = bcs_temp[bcid_n].temp_type
                            if ttype_n == temp_bc_heatflux:
                                bc_fluxn = bcs_temp[bcid_n].heat_flux
                                t_dd = t_d + bc_fluxn*dy
                            elif ttype_n == temp_bc_constant:
                                t_dd = fp(2.0)*bc_tn-t_d
                            rho = fp(0.5) * (dens[i,j,k] + dens[i,j+1,k])
                        elif bcs[bcid[i,j+1,k,fid_n]].type == bc_outlet:
                            bcid_n = bcid[i,j+1,k,fid_n]
                            bc_tn = bcs_temp[bcid_n].t
                            ttype_n = bcs_temp[bcid_n].temp_type
                            if ttype_n == temp_bc_heatflux:
                                bc_fluxn = bcs_temp[bcid_n].heat_flux
                                t_dd = t_d + bc_fluxn*dy
                            elif ttype_n == temp_bc_constant:
                                t_dd = fp(2.0)*bc_tn-t_d
                            rho = fp(0.5) * (dens[i,j,k] + dens[i,j+1,k])
                        else:
                            t_dd = t[i,k+2,j]
                            rho = fp(0.5) * (dens[i,j,k] + dens[i,j+1,k])
                        
                    f_d = rho*vf_d*area_y   # flux downstream
                    current_cell_boundary = False
                    
                    # if the south face of the current cell is boundary
                    if bc_s == bc_inlet:
                        #logging.debug("[South] Current cell is boundary")
                        if ttype_s == temp_bc_heatflux:
                            bc_fluxs = bcs_temp[bcid_s].heat_flux
                            t_u  = t_p + bc_fluxs*dy
                        elif ttype_s == temp_bc_constant:
                            bc_ts = bcs_temp[bcid_s].t
                            t_u = fp(2.0)*bc_ts-t_p 
                        t_uu = fp(2.0)*t_u-t_p
                        rho = dens[i,j,k]
                        current_cell_boundary = True
                    elif bc_s == bc_wall:
                        #logging.debug("[South] Current cell is boundary")
                        if ttype_s == temp_bc_heatflux:
                            bc_fluxs = bcs_temp[bcid_s].heat_flux
                            t_u  = t_p + bc_fluxs*dy
                        elif ttype_s == temp_bc_constant:
                            bc_ts = bcs_temp[bcid_s].t
                            t_u = fp(2.0)*bc_ts-t_p 
                        t_uu = fp(2.0)*t_u-t_p
                        rho = dens[i,j,k]
                        current_cell_boundary = True
                    elif bc_s == bc_outlet:
                        #logging.debug("[South] Current cell is boundary")
                        if ttype_s == temp_bc_heatflux:
                            bc_fluxs = bcs_temp[bcid_s].heat_flux
                            t_u  = t_p + bc_fluxs*dy
                        elif ttype_s == temp_bc_constant:
                            bc_ts = bcs_temp[bcid_s].t
                            t_u = fp(2.0)*bc_ts-t_p 
                        t_uu = fp(2.0)*t_u-t_p
                        rho = dens[i,j,k]
                        current_cell_boundary = True
                    else:
                        t_u = t[i,k-1,j]
                        rho = fp(0.5) * (dens[i,j-1,k] + dens[i,j,k])
                    
                    # if the south face of the south cell is boundary
                    if current_cell_boundary == False:
                        if bcs[bcid[i,j-1,k,fid_s]].type == bc_inlet:
                            bcid_s = bcid[i,j-1,k,fid_s]
                            bc_ts = bcs_temp[bcid_s].t
                            ttype_s = bcs_temp[bcid_s].temp_type
                            if ttype_s == temp_bc_heatflux:
                                bc_fluxs = bcs_temp[bcid_s].heat_flux
                                t_uu = t_u + bc_fluxs*dy
                            elif ttype_s == temp_bc_constant:
                                t_uu = fp(2.0)*bc_ts-t_u
                            rho = fp(0.5) * (dens[i,j,k] + dens[i,j-1,k])
                        elif bcs[bcid[i,j-1,k,fid_s]].type == bc_wall:
                            bcid_s = bcid[i,j-1,k,fid_s]
                            bc_ts = bcs_temp[bcid_s].t
                            ttype_s = bcs_temp[bcid_s].temp_type
                            if ttype_s == temp_bc_heatflux:
                                bc_fluxs = bcs_temp[bcid_s].heat_flux
                                t_uu = t_u + bc_fluxs*dy
                            elif ttype_s == temp_bc_constant:
                                t_uu = fp(2.0)*bc_ts-t_u
                            rho = fp(0.5) * (dens[i,j,k] + dens[i,j-1,k])
                        elif bcs[bcid[i,j-1,k,fid_s]].type == bc_outlet:
                            bcid_s = bcid[i,j-1,k,fid_s]
                            bc_ts = bcs_temp[bcid_s].t
                            ttype_s = bcs_temp[bcid_s].temp_type
                            if ttype_s == temp_bc_heatflux:
                                bc_fluxs = bcs_temp[bcid_s].heat_flux
                                t_uu = t_u + bc_fluxs*dy
                            elif ttype_s == temp_bc_constant:
                                t_uu = fp(2.0)*bc_ts-t_u
                            rho = fp(0.5) * (dens[i,j,k] + dens[i,j-1,k])
                        else:
                            t_uu = t[i,j-2,k]
                            rho = fp(0.5) * (dens[i,j-1,k] + dens[i,j,k])
                    
                    if ncy == 1:
                        t_uu = t_p; t_u = t_p; t_d = t_p; t_dd = t_p
                        rho  = dens[i,j,k]
                        
                    f_u = rho*vf_u*area_y   # flux upstream
                    
                    au = fp(0.); ad = fp(0.)
                    if f_u > fp(0.):
                        au = fp(1.)
                    if f_d > fp(0.):
                        ad = fp(1.)
                    
                    rdm = (t_dd - t_d) / (t_d - t_p + fp(1.e-12))    # rdm : r_downstream^minus
                    rdp = (t_p - t_u)  / (t_d - t_p + fp(1.e-12))     # rdp : r_downstream^plus
                    rum = (t_d - t_p)  / (t_p - t_u + fp(1.e-12))
                    rup = (t_u - t_uu) / (t_p - t_u + fp(1.e-12))
                    
                    if(t_dd==t_d) : rdm = fp(0.)
                    if(t_p==t_u) : rdp = fp(0.)
                    if(t_d==t_p) : rum = fp(0.)
                    if(t_u==t_uu) : rup = fp(0.)
                    
                    #logging.debug("t_p = " + str(t_p))
                    #logging.debug("t_d = " + str(t_d))
                    #logging.debug("rdp = " + str(rdp))
                    #logging.debug("rum = " + str(rum))
                    #logging.debug("rup = " + str(rup))
                    
                    src += (fp(0.5)*f_d*((fp(1.)-ad)*limiter(rdm,limiter_scheme)-ad*limiter(rdp,limiter_scheme))*(t_d-t_p) 
                                - fp(0.5)*f_u*((fp(1.)-au)*limiter(rum,limiter_scheme)-au*limiter(rup,limiter_scheme))*(t_p-t_u))
                    
                    ct[i,j,k,id_bsrc] += src
                    #logging.debug("SOU_src = " + str(src))
    return

