if True:
    import sys
    import numpy as np
    from numpy import sin, log, pi, angle, sqrt
    from scipy.special import sph_harm
    from scipy import special


    
    
def PAD_Momentum(l_m_list):

    x_momentum = np.linspace(-2.0 , 2.0, 40)
    y_momentum = np.linspace(-2.0 , 2.0, 40)
    z_momentum = np.linspace(-2.0 , 2.0, 40)
 
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
               
                out_going_wave = 0.0
                for l_m in l_m_list:
                    l, m  = l_m[0], l_m[1]
                    out_going_wave += sph_harm(m, l, phi, theta)

                pad_value_save[i,j,l] = out_going_wave


    pad_value_save = pad_value_save / pad_value_save.max()
    arr_reshaped = pad_value_save.reshape(pad_value_save.shape[0], -1)
    
    np.savetxt("SH.txt", arr_reshaped)

if __name__=="__main__":
    
    l_m_list = [[1,0]]
    PAD_Momentum(l_m_list)