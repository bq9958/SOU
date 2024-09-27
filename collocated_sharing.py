# -*- coding: utf-8 -*-
"""
Created on 04/05/2023
@author: 曾导SJTU
"""

import numpy as np
from precision_data import fp
from kernels import *
import logging
logger = logging.getLogger(__name__)

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
                    rho_e   = dens[i,j,k]
                else:
                    rho_e   = fp(0.5)*(dens[i,j,k]+dens[i+1,j,k])

                mul = con; mur = con
                ul    = uf[i+1,j,k]   # uniform velocity field, no interpolation needed
                ur    = uf[i+1,j,k]
                u_e = fp(0.5)*(ul+ur)
                aE    = a_nb(conv_scheme, area_x, idx, ul, ur, mul, mur, rho_e, -fp(1.0))

                # West Face
                if(i==0):
                    rho_w   = dens[i,j,k]
                else:
                    rho_w   = fp(0.5)*(dens[i,j,k]+dens[i-1,j,k])

                mul = con; mur = con
                ul    = uf[i,j,k]
                ur    = uf[i,j,k]
                u_w = fp(0.5)*(ul+ur)
                aW    = a_nb(conv_scheme, area_x, idx, ul, ur, mul, mur, rho_w, fp(1.0))

                # North Face
                if(j==ncy-1):
                    rho_n   = dens[i,j,k]
                else:
                    rho_n   = fp(0.5)*(dens[i,j,k]+dens[i,j+1,k])

                mul = con; mur = con
                ul  = vf[i,j+1,k]
                ur  = vf[i,j+1,k]
                u_n = fp(0.5)*(ul+ur)
                aN  = a_nb(conv_scheme, area_y, idy, ul, ur, mul, mur, rho_n, -fp(1.0))

                # South Face
                if(j==0):
                    rho_s   = dens[i,j,k]
                else:
                    rho_s   = fp(0.5)*(dens[i,j,k]+dens[i,j-1,k])

                mul = con; mur = con
                ul  = vf[i,j,k]
                ur  = vf[i,j,k]
                u_s = fp(0.5)*(ul+ur)
                aS  = a_nb(conv_scheme, area_y, idy, ul, ur, mul, mur, rho_s, fp(1.0))

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
                if conv_scheme == 3:
                    fe = u_e*rho_e
                    fw = u_w*rho_w
                    fn = u_n*rho_n
                    fs = u_s*rho_s
                    aP = (aP + aE + aW + aN + aS + aT + aB - sP*vol +  
                            (fe - fw) + (fn - fs))   # 在SOU中先不考虑瞬态项, to extend to 3D
                else:
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
    for k in range(ncz):
        for j in range(ncy):
            for i in range(ncx):
                if i==0:
                    ct[i,j,k,id_aP]   = fp(1.0)
                    ct[i,j,k,id_aE]   = fp(0.0)
                    ct[i,j,k,id_aW]   = fp(0.0)
                    ct[i,j,k,id_aN]   = fp(0.0)
                    ct[i,j,k,id_aS]   = fp(0.0)
                    ct[i,j,k,id_bsrc] = bc_tw
                elif i==ncx-1:
                    ct[i,j,k,id_aP]   = fp(1.0)
                    ct[i,j,k,id_aE]   = fp(0.0)
                    ct[i,j,k,id_aW]   = fp(0.0)
                    ct[i,j,k,id_aN]   = fp(0.0)
                    ct[i,j,k,id_aS]   = fp(0.0)
                    ct[i,j,k,id_bsrc] = bc_te

    return

