# -*- coding: utf-8 -*-
"""
Created on 04/05/2023
@author: 曾导SJTU
"""

import numpy as np
from precision_data import fp
from collocated_sharing import *
from solver import *
import logging
import log_config

def collocated_segregated(it, case, fluid, fluidboundary, post):
    
    # store old solution
    case.t0 = np.copy(case.t)

    dim     = case.dim
    
    ncx     = case.ncx    
    ncy     = case.ncy
    ncz     = case.ncz
        
    ncoef   = case.ncoef  
    dt      = case.dt
    
    x       = case.x 
    y       = case.y  
    z       = case.z 
    
    t       = case.t
    t0      = case.t0   
    uf      = case.uf     
    vf      = case.vf     
    wf      = case.wf     
    ct      = case.ct  
    r       = case.r   
    
    niter_t = case.niter_t
    relax_t = case.relax_t
     
    ieqn = case.ieqn
    eqn_conduction = case.eqn_conduction
    temp_tol = case.temp_tol
    res_t = case.res_t
    stop_sim = case.stop_sim
    
    # define local variables related to the fluid object
    spht     = fluid.spht
    con      = fluid.con 
    dens     = fluid.dens
    heat_src = fluid.heat_src
    
    if ieqn == eqn_conduction:
        
        # Evaluate temperature equation coefs
        conduction_coefs(case, dim, ncx, ncy, ncz, ncoef, dt, spht, con, heat_src, x, y, z, dens, t, uf, vf, wf, ct)
        if case.conv_scheme == 3:
            SOU_src(case, fluidboundary, dim, ncx, ncy, ncz, t, uf, vf, wf, dens, ct)
        #logging.debug("SOU_sec has been added")
        conduction_coef_bcs(case, fluidboundary, dim, ncx, ncy, ncz, ncoef, dt, con, x, y, z, dens, t, ct)
        
        initzero = False
        scalar_pj(case, post, dim, it, niter_t, relax_t, ncx, ncy, ncz, ncoef, ct, t, initzero, res_t) 
        # scalar_gs(case, post, dim, it, niter_t, relax_t, ncx, ncy, ncz, ncoef, ct, t, initzero, res_t)
        
        (case.l2_t, case.l2_max_t) = eqn_scalar_norm2(case, dim, it, ncx, ncy, ncz, t0, t, 'temp')

        if case.l2_t/case.l2_max_t < temp_tol:
            stop_sim = True
            print('')
            print('----------------------------')
            print('Final iter =', it)
            print('it, l2_t/l2_max_t', it, case.l2_t/case.l2_max_t)
            print('----------------------------')
        
        return stop_sim