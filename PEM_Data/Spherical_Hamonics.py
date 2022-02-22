if True:
    import sys
    import json 
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use('Agg')
    from matplotlib.colors import LogNorm
    from numpy import sin, log, pi, angle, sqrt
    import numpy as np
    import mpmath as mp
    from math import floor
    from scipy.special import sph_harm
    from scipy import special
 


def Cont_State_Calculator(l, k, grid):
    z = 2

    coulomb_fun = np.zeros(len(grid))
    for i, r in enumerate(grid):
        coulomb_fun[i] = mp.coulombf(l, -z/k, k*r)

    return coulomb_fun

def Index_Map(l_max, m_max):
    block_to_qn = {}
    qn_to_block = {}
    block = 0

    for m in range(0, m_max + 1):
            if m > 0:
                m_minus = -1*m
                for l in range(abs(m_minus), l_max + 1):
                    block_to_qn[block] = (l,m_minus)
                    qn_to_block[(l,m_minus)] = block
                    block += 1
                    
            for l in range(m, l_max + 1):
                block_to_qn[block] = (l,m)
                qn_to_block[(l,m)] = block
                block += 1
    return  block_to_qn, qn_to_block


def closest(lst, k): 
    return lst[min(range(len(lst)), key = lambda i: abs(float(lst[i])-k))] 

def PAD_Momentum():

    block_to_qn, qn_to_block = Index_Map(5, 0)
    grid = np.arange(0.5, 10 + 0.5, 0.5)
    x_momentum = np.linspace(-2.0 , 2.0, 20)
    y_momentum = np.linspace(-2.0 , 2.0, 20)
    z_momentum = np.linspace(-2.0 , 2.0, 20)
  

    pad_value = np.zeros((z_momentum.size,x_momentum.size))
    pad_value_save = np.zeros((z_momentum.size,y_momentum.size, x_momentum.size))
    
    for i, px in enumerate(x_momentum):
        print(round(px,3))
        for j, py in enumerate(y_momentum):

            for l, pz in enumerate(z_momentum):

                k = np.sqrt(px*px + py*py + pz*pz)
                if k == 0:
                    continue

                if px > 0 and py > 0:
                    phi = np.arctan(py/px)
                elif px > 0 and py < 0:
                    phi = np.arctan(py/px) + 2*pi
                elif px < 0 and py > 0:
                    phi = np.arctan(py/px) + pi
                elif px < 0 and py < 0:
                    phi = np.arctan(py/px) + pi
                elif px == 0 and py == 0:
                    phi = 0
                elif px == 0 and py > 0:
                    phi = pi / 2
                elif px == 0 and py < 0:
                    phi = 3*pi / 2
                elif py == 0 and px > 0:
                    phi = 0
                elif py == 0 and px < 0:
                    phi = pi

                theta = np.arccos(pz/k)
                
                
                theta, phi = np.meshgrid(theta, phi)
                out_going_wave = np.zeros(phi.shape, dtype=complex)
                
                for key in qn_to_block:
                    l, m  = key[0], key[1]
                    
                    CF_Proj = Cont_State_Calculator(l, k, grid)
                    CF_Psi = Cont_State_Calculator(1, k, grid)
                    z=2
                    phase = angle(special.gamma(l + 1 - 1j*z/k))
                    coef =  np.exp(-1.0j*phase)* np.power(1.0j,l)*np.sum(CF_Proj.conj()*CF_Psi)  
                    out_going_wave += coef*sph_harm(m, l, phi, theta)


                pad_value_save[i,j,l] = np.abs(out_going_wave)**2
            


    pad_value_save = pad_value_save / pad_value_save.max()
    arr_reshaped = pad_value_save.reshape(pad_value_save.shape[0], -1)
    
    np.savetxt("PEM.txt", arr_reshaped)
    
    exit()
    return pad_value, x_momentum, z_momentum


if __name__=="__main__":

    PAD_Momentum()

