from matplotlib import pyplot as plt
from os.path import expanduser

from utilities import get_s, MC_sweep, get_samples
from observables import *


def plot_equilibration():
    N = 50
    k = 6
    equib = 100_000
    Ts = np.linspace(0.01, 1.2*Tc, k)

    fig, ax = plt.subplots(2, k, figsize=(16, 8))
    fig.suptitle("$N={}$\n".format(N) + "{} sweeps".format(equib))

    for i, T in enumerate(Ts):
        print(i)
        s = get_s(N)
        ax[0, i].set_title("$T = {:.3f}$".format(T))
        ax[0, i].imshow(s)
        for _ in range(equib):
            MC_sweep(s, N, T)
        ax[1, i].imshow(s)

    plt.tight_layout()
    plt.savefig(expanduser("~") + "/Desktop/test2"  + ".png")
    plt.close(fig)

def plot_observables():
    observables = {
        "E": energy_density,
        "|M|": absolute_magnetization
    }

    equib = 100_000
    n = 100_000
    m = 1
    k = 20
    l = len(observables)

    Ns = [8, 16]
    Ts = np.linspace(1.5, 1.2*Tc, k)

    fig, ax = plt.subplots(l, figsize=(10, 10))

    if l==1: ax = [ax,]

    for i, N in enumerate(Ns):
        print(i)
        samples = np.empty((k, l))
        for j, T in enumerate(Ts):
            samples[j] = get_samples(N, T, n, m, equib, observables)
        for j, key in enumerate(observables):
            ax[j].plot(Ts, samples[:, j], ".", label="$N ={}$".format(N))
            ax[j].legend()
            ax[j].set_title = "$" + key + "$"
    
    plt.tight_layout()
    plt.savefig(expanduser("~") + "/Desktop/test"  + ".png")
    plt.close(fig)

def plot_heat_capacity():
    observables = {
        "energy_density": energy_density,
        "energy_density_squared": energy_density_squared
    }

    samples = {}
    equib = 1000
    n = 100
    m = 2
    k = 10
    l = len(observables)

    Ns = [8, 16]
    Ts = np.linspace(0.1, 1.2*Tc, k)

    fig, ax = plt.subplots(l)

    if l == 1: ax = [ax,]
    

    for i, N in enumerate(Ns):
        print(i)
        samples = np.empty((k, l))
        for j, T in enumerate(Ts):
            samples[j] = get_samples(N, T, n, m, equib, observables)

        heat_cap = samples[:, 1] - samples[:, 0]**2 
        ax[j].plot(Ts, heat_cap, ".", label="$N ={}$".format(N))

    plt.tight_layout()
    plt.savefig(expanduser("~") + "/Desktop/test3"  + ".png")
    plt.close(fig)
    

if __name__ == "__main__":
    # plot_equilibration()
    plot_observables()