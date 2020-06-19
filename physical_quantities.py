import numpy as np


# Critical temprature of the 2D ising model
Tc = 2 / (np.log(1 + np.sqrt(2)))

def absolute_magnetization(s, N, T):
    return np.abs(np.sum(s)) / N ** 2
    
def magnetization_squared(s, N, T):
    return (np.sum(s) / N**2)**2

def energy(s, N, T):
    sum_neigh = np.zeros_like(s)
    for j in range(2):
        for n in (-1, 1):
            sum_neigh += np.roll(s, n, axis=j)
    
    return -1/2 *np.sum(s * sum_neigh)

def energy_density(s, N, T):
    return energy(s, N, T) / N**2
    
def energy_density_squared(s, N, T):
    return energy_density(s, N, T)** 2

def fraction(s, N, T):
    energy_diff = 4*np.sum(-s[-1]*s[0])
    return np.exp(energy_diff/T)


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
    "energy_density": energy_density,
    "energy_density_squared": energy_density_squared,
    "absolute_magnetization": absolute_magnetization,
    "magnetization squared": magnetization_squared,
    "fraction": fraction
}

def heat_capacity(samples, T, N):
    return (samples["energy_density_squared"] - samples["energy_density"]**2) * np.outer(N, 1/T)**2

def susceptability(samples, T, N):
    return (samples["magnetization squared"] - samples["absolute_magnetization"]** 2) * np.outer(N ** 2, 1 / T)
    
def tension(samples, T, N):
    return - np.outer(1 / N, T) * np.log(samples["fraction"])
    
def reudced_temp(Ts):
    return (Tc - Ts) / Tc

func_of_obs = {
    "heat_capacity": heat_capacity,
    "susceptibility": susceptability,
    "tension": tension
}

units = {
    "energy_density":           "J",
    "energy_density_squared":   "J^2",
    "absolute_magnetization":   "Am^{-1}",
    "magnetization squared":    "A^2m^{-2}",
    "fraction":                 "1",
    "heat_capacity":            "1",
    "susceptibility":           "Am^{-1}J^{-1}",
    "tension":                  "J"
}

latex = {
    "energy_density": "E",
    "energy_density_squared": "E^2",
    "absolute_magnetization": "|M|",
    "magnetization squared": "M^2",
    "fraction": "Z_k/Z_t",
    "heat_capacity": "C_M",
    "susceptibility": "\chi",
    "tension" : "\\tau"
}