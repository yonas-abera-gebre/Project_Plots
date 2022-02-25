import scipy as sci
import scipy.special as sp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors



def Plot_SH():
    # Coordinate arrays for the graphical representation
    x = np.linspace(-np.pi, np.pi, 100)
    y = np.linspace(-np.pi/2, np.pi/2, 50)
    X, Y = np.meshgrid(x, y)

    # Spherical coordinate arrays derived from x, y
    # Necessary conversions to get Mollweide right
    phi = x.copy()    # physical copy
    phi[x < 0] = 2 * np.pi + x[x<0]
    theta = np.pi/2 - y
    PHI, THETA = np.meshgrid(phi, theta)

    l = 2
    m = 2
    SH_SP = sp.sph_harm(m, l, PHI, THETA)

    l = 2
    m = -2
    SH_SP += sp.sph_harm(m, l, PHI, THETA)
    
    
    l = 2
    m = 0
    SH_SP += sp.sph_harm(m, l, PHI, THETA) * 2
    
    # l = 4
    # m = 2
    # SH_SP += sp.sph_harm(m, l, PHI, THETA)*0.2
    
    # l = 4
    # m = -2
    # SH_SP += sp.sph_harm(m, l, PHI, THETA)*0.2
     
    # l = 4
    # m = 0
    # SH_SP += sp.sph_harm(m, l, PHI, THETA) * 0.5
    
    SH_SP = np.abs(SH_SP)
    
    matplotlib.rc('text', usetex=True)
    matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
    
    xlabels = ['$210^\circ$', '$240^\circ$','$270^\circ$','$300^\circ$','$330^\circ$',
            '$0^\circ$', '$30^\circ$', '$60^\circ$', '$90^\circ$','$120^\circ$', '$150^\circ$']

    ylabels = ['$165^\circ$', '$150^\circ$', '$135^\circ$', '$120^\circ$', 
            '$105^\circ$', '$90^\circ$', '$75^\circ$', '$60^\circ$',
            '$45^\circ$','$30^\circ$','$15^\circ$']

    fig, ax = plt.subplots(subplot_kw=dict(projection='mollweide'), figsize=(10,8))
    im = ax.pcolormesh(X, Y , SH_SP/np.max(SH_SP), cmap='viridis')
    ax.set_xticklabels(xlabels, fontsize=14, color='r')
    ax.set_yticklabels(ylabels, fontsize=14)
    # ax.set_title('real$(Y^2_ 4)$', fontsize=20)
    ax.set_xlabel(r'$\boldsymbol \phi$', fontsize=20)
    ax.set_ylabel(r'$\boldsymbol{\theta}$', fontsize=20)
    ax.grid()
    # fig.colorbar(im, orientation='horizontal')
    # plt.savefig("Two_Photon_Bonding.png")
    plt.show()

if __name__=="__main__":

    Plot_SH()

