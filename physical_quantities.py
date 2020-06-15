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

    return - 1 / 2 * np.sum(s * sum_neigh) / N ** 2
    
def energy_density_squared(s, N):
    return energy_density(s, N)** 2


observables = {
    "E": energy_density,
    "E2": energy_density_squared,
    "abs_M": absolute_magnetization,
    "M2": magnetization_squared,
}
