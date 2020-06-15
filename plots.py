import numpy as np
from matplotlib import pyplot as plt
from os import path, mkdir

from utilities import get_s, MC_sweep, get_samples, read_samples
from physical_quantities import observables, Tc


def plot_equilibration():
    N = 200
    k = 6
    equib = 100
    Ts = np.linspace(0.01, 1.2*Tc, k)

    fig, ax = plt.subplots(2, k, figsize=(24, 12))
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
    plt.savefig("figs/equilibration.png")
    plt.close(fig)


def plot_observables(sub_dir=""):
    Ns, Ts, samples = read_samples(list(observables), sub_dir=sub_dir)
    fig_path = "figs/" + sub_dir
    if not path.isdir(fig_path):
        mkdir(fig_path)

    for quantity in list(observables):
        fig, ax = plt.subplots(figsize=(10, 8))
        fig.suptitle = quantity

        for i, N in enumerate(Ns):
            ax.plot(Ts, samples[quantity][i], ".", label="$N ={}$".format(N))
            ax.legend()
        
        plt.tight_layout()
        plt.savefig(fig_path + quantity + ".png")
        plt.close(fig)

