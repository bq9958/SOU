# -*- coding: utf-8 -*-
"""
Created on 04/05/2023
@author: 曾导SJTU
"""

import numpy as np
from precision_data import fp

def eqn_scalar_norm2(case, dim, it_nl, ncx, ncy, ncz, u0, u, var):

    # define local variables related to the case object
    if var == "temp":
        l2_u = case.l2_t
        l2_max_u = case.l2_max_t

    # Compute the L2 norm
    
    # Equivalen form
    # for k in range(ncz):
    #     for j in range(ncy):
    #         for i in range(ncx):
    #             res_u = u[i, j, k] - u0[i, j, k]
    #             l2_u += res_u**2
    # ncell = ncx*ncy*ncz
    # l2_u = math.sqrt(l2_u/fp(ncell))
    
    l2_u = np.sqrt(np.mean((u - u0)**2))

    # Compute the maximum L2 norm observed so far
    l2_max_u = max(l2_u, l2_max_u)

    return (l2_u, l2_max_u)

def scalar_pj(case, post, dim, it_nl, niter, relax, ncx, ncy, ncz, ncoef, ct, t, initzero, res):
    
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
    
    nsteps = case.nsteps
    
    # define local variables related to the post object
    res_freq = post.res_freq
    linsol_fname = post.linsol_fname
    
    if initzero:
        t.fill(fp(0.0))

    # Allocate t0 array and initialize it
    t0 = np.zeros((ncx, ncy, ncz), dtype=fp)

    max_norm = -1.e20

    # Begin linear solver iterations
    for it in range(1, niter+1):
        
        # Copy current to old
        t0 = t.copy()

        # Initialize norm
        norm2 = 0.0

        for k in range(ncz):
            for j in range(ncy):
                for i in range(ncx):
                    
                    if i == 0:
                        tw = 0.0
                    else:
                        tw = t0[i-1,j,k]
                    if i == ncx-1:
                        te = 0.0
                    else:
                        te = t0[i+1,j,k]
                    
                    if j == 0:
                        ts = 0.0
                    else:
                        ts = t0[i,j-1,k]
                    if j == ncy-1:
                        tn = 0.0
                    else:
                        tn = t0[i,j+1,k]
                    
                    if dim == 3:
                        if k == 0:
                            tb = 0.0
                        else:
                            tb = t0[i,j,k-1]
                        if k == ncz-1:
                            tt = 0.0
                        else:
                            tt = t0[i,j,k+1]

                    t_new = (                 \
                        - ct[i,j,k,id_aE]*te  \
                        - ct[i,j,k,id_aW]*tw  \
                        - ct[i,j,k,id_aN]*tn  \
                        - ct[i,j,k,id_aS]*ts  \
                        + ct[i,j,k,id_bsrc]   \
                            )

                    if dim==3:
                        t_new = t_new - ct[i,j,k,id_aT]*tt - ct[i,j,k,id_aB]*tb

                    t_new = t_new / ct[i,j,k,id_aP]

                    du = relax * (t_new - t0[i,j,k])
                    t[i,j,k] = t0[i,j,k] + du
                    norm2 = norm2 + du*du

        ncell = ncx*ncy*ncz
        norm2 = np.sqrt(norm2 / fp(ncell))

        if it == 1:
            initial_norm = norm2+1.e-20

        max_norm = max(norm2, max_norm) + 1.e-20

        rel_norm = norm2 / max_norm

        if rel_norm < res or it == niter:
            case.total_linsol_iters = case.total_linsol_iters + it
            if it_nl % res_freq == 0 or it_nl == 1 or it_nl == nsteps:
                print('it_nl, it, tot_it, norm2, init, max, rel ', it_nl, it, case.total_linsol_iters, norm2, initial_norm, max_norm, rel_norm)
                with open(linsol_fname, 'a') as linsol_fid:
                    print(it_nl, it, case.total_linsol_iters, norm2, initial_norm, max_norm, rel_norm, file=linsol_fid)
            break

    return

def scalar_gs(case, post, dim, it_nl, niter, relax, ncx, ncy, ncz, ncoef, ct, t, initzero, res):
    
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
    
    nsteps = case.nsteps
    
    # define local variables related to the post object
    res_freq = post.res_freq
    linsol_fname = post.linsol_fname
    
    if initzero:
        t.fill(fp(0.0))

    max_norm = -1.e20

    # Begin linear solver iterations
    for it in range(1, niter+1):
        
        # Initialize norm
        norm2 = 0.0

        for k in range(ncz):
            for j in range(ncy):
                for i in range(ncx):
                    
                    if i == 0:
                        tw = 0.0
                    else:
                        tw = t[i-1,j,k]
                    if i == ncx-1:
                        te = 0.0
                    else:
                        te = t[i+1,j,k]
                    
                    if j == 0:
                        ts = 0.0
                    else:
                        ts = t[i,j-1,k]
                    if j == ncy-1:
                        tn = 0.0
                    else:
                        tn = t[i,j+1,k]
                    
                    if dim == 3:
                        if k == 0:
                            tb = 0.0
                        else:
                            tb = t[i,j,k-1]
                        if k == ncz-1:
                            tt = 0.0
                        else:
                            tt = t[i,j,k+1]

                    t_new = (                 \
                        - ct[i,j,k,id_aE]*te  \
                        - ct[i,j,k,id_aW]*tw  \
                        - ct[i,j,k,id_aN]*tn  \
                        - ct[i,j,k,id_aS]*ts  \
                        + ct[i,j,k,id_bsrc]   \
                            )

                    if dim==3:
                        t_new = t_new - ct[i,j,k,id_aT]*tt - ct[i,j,k,id_aB]*tb

                    t_new = t_new / ct[i,j,k,id_aP]

                    du = relax * (t_new - t[i,j,k])
                    t[i,j,k] = t[i,j,k] + du
                    norm2 = norm2 + du*du

        ncell = ncx*ncy*ncz
        norm2 = np.sqrt(norm2 / fp(ncell))

        if it == 1:
            initial_norm = norm2+1.e-20

        max_norm = max(norm2, max_norm) + 1.e-20

        rel_norm = norm2 / max_norm

        if rel_norm < res or it == niter:
            case.total_linsol_iters = case.total_linsol_iters + it
            if it_nl % res_freq == 0 or it_nl == 1 or it_nl == nsteps:
                print('it_nl, it, tot_it, norm2, init, max, rel ', it_nl, it, case.total_linsol_iters, norm2, initial_norm, max_norm, rel_norm)
                with open(linsol_fname, 'a') as linsol_fid:
                    print(it_nl, it, case.total_linsol_iters, norm2, initial_norm, max_norm, rel_norm, file=linsol_fid)
            break

    return

#    def scalar_adi ...
