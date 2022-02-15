import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('Agg')
matplotlib.use( 'tkagg' )
from numpy import sin, cos, pi, exp, sqrt, dot
from scipy import integrate
import scipy.special as special
import Module as Mod 
import mpmath as mp
from scipy.signal import find_peaks
from math import floor, ceil
import plotly.graph_objects as go
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
    return round(N)

def Calculate_NPA2(k, quiver_energy, Ip, omega):
    N = (pow(k,2)/2 + quiver_energy + Ip)/omega
    return ceil(N)

def GBF_Zeros(k_array, quiver_energy, quiver_radius, Ip, omega):
    
    BF_array = np.zeros(len(k_array))
    
    for i, k in enumerate(k_array):
        N = Calculate_NPA(k, quiver_energy, Ip, omega)
        BF_array[i] = k*pow(GBF(N, quiver_radius*k, quiver_energy/(2*omega)),2)*pow(quiver_energy - N*omega, 2)
    
    return BF_array

def GBF_Zeros2(k_array, quiver_energy, quiver_radius, Ip, omega):
    
    BF_array = np.zeros(len(k_array))
    
    for i, k in enumerate(k_array):
        N = Calculate_NPA2(k, quiver_energy, Ip, omega)
        BF_array[i] = k*pow(GBF(N, quiver_radius*0, quiver_energy/(2*omega)),2)*pow(quiver_energy - N*omega, 2)
    
    return BF_array

    
if __name__=="__main__":
    

    X, Y, Z = np.mgrid[-2:2:80j, -2:2:80j, -2:2:80j]
    
    x_momentum = np.linspace(-2.0 , 2.0, 80)
    y_momentum = np.linspace(-2.0 , 2.0, 80)
    z_momentum = np.linspace(-2.0 , 2.0, 80)
    arr = np.zeros((z_momentum.size,y_momentum.size, x_momentum.size))
    
    file_name = sys.argv[1]
    loaded_arr = np.loadtxt(file_name)
    load_original_arr = loaded_arr.reshape(loaded_arr.shape[0], loaded_arr.shape[1] // arr.shape[2], arr.shape[2])
    

    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=load_original_arr.flatten(),
        isomin=0.0,
        isomax=1.0,
        opacity=0.1, # needs to be small to see through all surfaces
        surface_count=17, # needs to be a large number for good volume rendering
        ))
    fig.show()
   
    exit()
    
    intensity = 3.51e16
    omega = 0.75
    Ip = 1.0
    k_array = np.arange(0, 4, 0.05)
    e_array = np.power(k_array,2)/2.0

    quiver_radius, quiver_energy = Calculate_Parameters(intensity, omega)
    
    EN_values = [N*omega - quiver_energy - Ip for N in range(10)]
    
    KN_values = [sqrt(2*en) for en in EN_values]
    

    BF_array = GBF_Zeros2(k_array, quiver_energy, quiver_radius, Ip, omega)
    # BF_array /= np.max(BF_array)
    
    
    plt.plot(k_array, BF_array)
    
    for kn in KN_values:
        plt.axvline(x=kn, color='r')
    
    # for x in np.arange(0, 3, 0.5):
    #     plt.axvline(x=x, color='b')

    
    plt.xlim(0, 4.0)
    # plt.ylim(1e-2, 1)
    # plt.yscale("log")
    plt.show()
    plt.savefig("BF_K2.png")
        