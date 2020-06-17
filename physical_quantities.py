import numpy as np


# Critical temprature of the 2D ising model
Tc = 2 / (np.log(1 + np.sqrt(2)))

def absolute_magnetization(s, N):
    return np.abs(np.sum(s)) / N ** 2
    
def magnetization_squared(s, N):
    return (np.sum(s) / N**2)**2

def energy_density(s, N):
    sum_neigh = np.zeros_like(s)
    for j in range(2):
        for n in (-1, 1):
            sum_neigh += np.roll(s, n, axis=j)

    return - np.sum(s * sum_neigh) / N ** 2
    
def energy_density_squared(s, N):
    return energy_density(s, N)** 2

def energy_MJ(s, N, border=(1, 1)):
    """ gives energy as in the Mon-Jasnow scheme """
    
    row = np.ones((1, N), dtype=int)
    s_w_border = np.concatenate([border[0]*row, s, border[1]*row])
    sum_neigh = np.zeros_like(s)

    for j in range(2):
        for n in (-1, 1):
            sum_neigh += np.roll(s_w_border, n, axis=j)[1:-1]

    return - 1/2*np.sum(s * sum_neigh)
    
def ms(s, N):
    return np.sum(s[-1])


observables = {
    "E": energy_density,
    "E^2": energy_density_squared,
    "|M|": absolute_magnetization,
    "M^2": magnetization_squared,
}

def heat_capacity(samples, T, N):
    return (samples["E^2"] - samples["E"]**2) * np.outer(N, 1/T)**2

def susceptability(samples, T, N):
    return (samples["M^2"] - samples["|M|"]**2) * np.outer(N**2, 1/T)

func_of_obs = {
    "C_M": heat_capacity,
    "\chi": susceptability
}

units = {
    "E": "J",
    "E^2": "J^2",
    "|M|": "Am^{-1}",
    "M^2": "A^2m^{-2}",
    "C_M": "1",
    "\chi": "Am^{-1}J^{-1}"
}

name = {
    "E": "energy_density",
    "E^2": "energy_density_squared",
    "|M|": "absolute_magnetization",
    "M^2": "magnetization squared",
    "C_M": "heat_capacity",
    "\chi": "susceptibility"  
}