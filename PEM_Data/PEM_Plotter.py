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
import plotly.graph_objects as go
from plotly.subplots import make_subplots
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


def File_Name_To_Data(file_name, X_lim = (-2,2,80), Y_lim = (-2,2,80), Z_lim = (-2,2,80)):
    
    x_momentum = np.linspace(X_lim[0],X_lim[1], X_lim[2])
    y_momentum = np.linspace(Y_lim[0],Y_lim[1], Y_lim[2])
    z_momentum = np.linspace(Z_lim[0],Z_lim[1], Z_lim[2])
    arr = np.zeros((z_momentum.size,y_momentum.size, x_momentum.size))
    
    
    loaded_arr = np.loadtxt(file_name)
    PEM_data = loaded_arr.reshape(loaded_arr.shape[0], loaded_arr.shape[1] // arr.shape[2], arr.shape[2])
    
    return  PEM_data

def PEM_Data_Plotter(PEM_data, X_lim = (-2,2,80), Y_lim = (-2,2,80), Z_lim = (-2,2,80)):
    
    X, Y, Z = np.mgrid[X_lim[0]:X_lim[1]:X_lim[2]*1.0j, Y_lim[0]:Y_lim[1]:Y_lim[2]*1.0j, Z_lim[0]:Z_lim[1]:Z_lim[2]*1.0j]


    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=PEM_data.flatten(),
        isomin=0.0,
        isomax=1.0,
        opacity=0.1, # needs to be small to see through all surfaces
        surface_count=17, # needs to be a large number for good volume rendering
        ))
    fig.show()

def PEM_Data_Plotter_Multiple(PEM_data_list, X_lim = (-2,2,80), Y_lim = (-2,2,80), Z_lim = (-2,2,80)):
    fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'volume'}, {'type': 'volume'}]])

    X, Y, Z = np.mgrid[X_lim[0]:X_lim[1]:X_lim[2]*1.0j, Y_lim[0]:Y_lim[1]:Y_lim[2]*1.0j, Z_lim[0]:Z_lim[1]:Z_lim[2]*1.0j]
    
    fig.add_trace(go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=PEM_data_list[0].flatten(),
        isomin=0.0,
        isomax=1.0,
        opacity=0.1, # needs to be small to see through all surfaces
        surface_count=17, # needs to be a large number for good volume rendering
        ), row=1, col=1)
    fig.add_trace(go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=PEM_data_list[0].flatten(),
        isomin=0.0,
        isomax=1.0,
        opacity=0.1, # needs to be small to see through all surfaces
        surface_count=17, # needs to be a large number for good volume rendering
        ), row=1, col=2)
    
    # fig.update_traces(x=X.flatten(), y=Y.flatten(), z=Z.flatten(), value=values.flatten(),isomin=0.15, isomax=0.9, opacity=0.1, surface_count=15)
    fig.show()
    
if __name__=="__main__":
    
    if len(sys.argv) == 2:
        file_name = sys.argv[1]
        PEM_data = File_Name_To_Data(file_name)
        PEM_Data_Plotter(PEM_data)
    if len(sys.argv) == 3:
        file_name_1 = sys.argv[1]
        PEM_data_1 = File_Name_To_Data(file_name_1)
        file_name_2 = sys.argv[2]
        PEM_data_2 = File_Name_To_Data(file_name_2)
        
        PEM_data_list = [PEM_data_1, PEM_data_2]
    
        PEM_Data_Plotter_Multiple(PEM_data_list, X_lim = (-2,2,80), Y_lim = (-2,2,80), Z_lim = (-2,2,80))
    
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
        