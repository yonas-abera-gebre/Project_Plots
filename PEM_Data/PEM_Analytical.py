import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('Agg')
matplotlib.use( 'tkagg' )
from numpy import sin, cos, pi, exp, sqrt, dot
from scipy import integrate
import scipy.special as special
import mpmath as mp
from scipy.signal import find_peaks
from math import floor, ceil
# import plotly.graph_objects as go
# from plotly.subplots import make_subplots
import sys

def GBF(n,a,b):

    def func_for_integration(phi):
        return_val = 1.0j*(a*sin(phi) + b*sin(2*phi) - n*phi)
        return exp(return_val)

    integral = integrate.quadrature(func_for_integration, -1.0*pi, pi)

    if abs(integral[0].imag) > 1e-5:
        print("imag part large")
        exit()
    else:
        return 1/(2*pi)*integral[0].real

def GBF_2(n,a,b):

    return_val = 0.0j
    for k in range(-50,50,1):
        return_val += special.jv(n-int(2*k), a)*special.jv(k,b)
    print(return_val)
    return return_val

def GBF_3(n,a,b):
    
    return_val = 0.0j
    for k in range(-50,50,1):
        return_val += mp.besselj(n-int(2*k), a)*special.jv(k,b)
    return return_val.real

def Calculate_Parameters(intensity, omega):
    intensity = intensity/3.51e16
    quiver_radius = sqrt(intensity)/pow(omega,2)
    quiver_energy = intensity/(4*pow(omega,2))

    return quiver_radius, quiver_energy
    
def Calculate_NPA(k, quiver_energy, Ip, omega):
    N = (pow(k,2)/2 + quiver_energy + Ip)/omega
    return ceil(N)

def GBF_Zeros(k_array, quiver_energy, quiver_radius, Ip, omega):
    
    BF_array = np.zeros(len(k_array))
    
    for i, k in enumerate(k_array):
        N = Calculate_NPA2(k, quiver_energy, Ip, omega)
        BF_array[i] = k*pow(GBF(N, quiver_radius*0, quiver_energy/(2*omega)),2)*pow(quiver_energy - N*omega, 2)
    
    return BF_array

def PAD_2D(quiver_energy, quiver_radius, omega, polarization):


    x_momentum = np.linspace(-2.0 , 2.0, 80)
    y_momentum = np.linspace(-2.0 , 2.0, 80)
    z_momentum = np.linspace(-2.0 , 2.0, 80)

    
    pad_value = np.zeros((z_momentum.size,x_momentum.size))
    pad_value_save = np.zeros((z_momentum.size,y_momentum.size, x_momentum.size))
     
    for i, kx in enumerate(x_momentum):
        print(round(kx,3))
        for j, kz in enumerate(z_momentum):
            pad_value_temp = 0.0
            for l, ky in enumerate(y_momentum):
                
                k = np.sqrt(kx*kx + ky*ky + kz*kz)
                N = Calculate_NPA(k, quiver_energy, Ip, omega)
                k_vector= np.array([kx,ky,kz], dtype=float)
             
                mole_orientation = np.array([0,0,1])
                quiver_vector = quiver_radius*polarization
                
                pad_value_temp += k*pow(GBF(N, np.dot(quiver_vector, k_vector), quiver_energy/(2*omega)),2)*pow(quiver_energy - N*omega ,2)*pow(np.cos(np.dot(k_vector, mole_orientation)), 2)

                pad_value_save[i,j,l] = k*pow(GBF(N, np.dot(quiver_vector, k_vector), quiver_energy/(2*omega)),2)*pow(quiver_energy - N*omega ,2)*pow(np.cos(np.dot(k_vector, mole_orientation)), 2)
                
            pad_value[j, i] += pad_value_temp
    
    pad_value_save = pad_value_save / pad_value_save.max() 
    arr_reshaped = pad_value_save.reshape(pad_value_save.shape[0], -1)
    np.savetxt("Ana_X_B_One_Photon.txt", arr_reshaped)
    
    exit()
    
    pad_value = pad_value / pad_value.max()
    pos = plt.imshow(pad_value, extent=[-2.0, 2.0, -2.0, 2.0])#, norm=LogNorm(vmin=1e-3, vmax=1))
    plt.colorbar(pos)
    plt.savefig("FIG.png")
    
if __name__=="__main__":
    
    polarization = np.array([1,0,0])
    ellipticity = 0.
    ellipticity_vector = np.array([0,1,0])
    intensity = 1e14
    omega = 1.5
    Ip = 1.0
    quiver_radius, quiver_energy = Calculate_Parameters(intensity, omega)
    PAD_2D(quiver_energy, quiver_radius, omega, polarization)
    
