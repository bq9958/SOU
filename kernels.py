# -*- coding: utf-8 -*-
"""
Created on 04/05/2023
@author: 曾导SJTU
"""

import numpy as np
from precision_data import fp
import logging
logger = logging.getLogger(__name__)

def cal_area_vol(dx, dy, dz):
    
    area_x = dy*dz
    area_y = dx*dz
    area_z = dx*dy
    vol    = dx*area_x
    
    return area_x, area_y, area_z, vol

def a_pec_pow(pec):
    # Incoming variable
    # pec: float

    ap = fp(1.0) - fp(0.1) * abs(pec)
    ap = max(fp(0.0), ap**5)

    return ap

def a_nb(conv_scheme, area, idx, f, gl, gr, sign_f):

    d = fp(2.0)*gl*gr/(gl+gr+fp(1.e-12))*idx

    # print("conv_scheme ", conv_scheme)
    
    if conv_scheme == 0:
        # Upwind
        # Pantakar Table 5.2
        a = area * ( d + max(fp(0.0),sign_f*f) )
        # print('a ', a)
    elif conv_scheme == 1:
        # Central Difference
        # Pantakar Table 5.2中的CD scheme和红宝书Eq. (11.25)等价，课上推导
        a = area * ( d*(fp(1.0)-fp(0.5)*abs(f/d)) + max(fp(0.0),sign_f*f) )
        # print('a ', a)
    elif conv_scheme == 2:
        # Power-law
        # Pantakar Table 5.2
        a = area * ( d*a_pec_pow(abs(f/d)) + max(fp(0.0),sign_f*f) )
    elif conv_scheme == 3:
        # SOU
        a = area * ( d + max(fp(0.0),sign_f*f) )  # same aw, ae, as, an as Upwind
    
    return a

def limiter(limiter,r):
    if limiter == 0:   #vanLeer limiter
        if r == fp(-1.0):
            return r   # why output r ? 
        else:
            return (r+abs(r))/(1+r)
    
def SOU_source(limiter, fw, fe, fn, fs):
    if fw > 0:
        alpha_w = 1
    else:
        alpha_w = 0
    if fe > 0:
        alpha_e = 1
    else:
        alpha_e = 0
    if fs > 0:
        alpha_s = 1
    else:
        alpha_s = 0
    if fn > 0:
        alpha_n = 1
    else:
        alpha_n = 0     

                        
    
        
    