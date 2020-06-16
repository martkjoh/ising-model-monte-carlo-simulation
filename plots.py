import numpy as np
from matplotlib import pyplot as plt
from os import path, mkdir

from utilities import get_s, MC_sweep, get_samples, read_samples
from physical_quantities import observables, func_of_obs, Tc, units, name


font = {
    'family': 'serif',
    'weight': 'normal',
    'size': 30}
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rc("lines", linewidth=1, markersize=3)
styles = ["--x", "--o", "--+", "--d"]


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


def plot_observables(sub_dir):
    Ns, Ts, samples = read_samples(list(observables), sub_dir=sub_dir)
    fig_path = "figs/" + sub_dir
    if not path.isdir(fig_path):
        mkdir(fig_path)

    for quantity in list(observables):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_ylabel("$" + quantity + "/[\\mathrm{" + units[quantity] + "}]$")
        ax.set_xlabel("$T / [J]$")

        for i, N in enumerate(Ns):
            ax.plot(Ts, samples[quantity][i] / 2, styles[i], label="$N={}$".format(int(N)))
        ax.plot([Tc, Tc], ax.get_ylim(), "k--", label="$T_c$")
        ax.grid(True)
        ax.legend()        
        
        plt.tight_layout()
        plt.savefig(fig_path + name[quantity] + ".png", dpi=300)
        plt.close(fig)

def plot_funcs(sub_dir=""):
    Ns, Ts, samples = read_samples(list(observables), sub_dir=sub_dir)
    fig_path = "figs/" + sub_dir
    if not path.isdir(fig_path):
        mkdir(fig_path)

    for quantity in list(func_of_obs):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_ylabel("$" + quantity + "/[\\mathrm{" + units[quantity] + "}]$")
        ax.set_xlabel("$T / [J]$")

        y = func_of_obs[quantity](samples, Ts, Ns)
        for i, N in enumerate(Ns):
            ax.plot(Ts, y[i] / 2, styles[i], label="$N={}$".format(int(N)))
        ax.plot([Tc, Tc], ax.get_ylim(), "k--", label="$T_c$")
        ax.grid(True)
        ax.legend()

        plt.tight_layout()
        plt.savefig(fig_path + name[quantity] + ".png", dpi=300)
        plt.close(fig)

def mag_plot(sub_dir):
    Ns, Ts, samples = read_samples(list(observables), sub_dir=sub_dir)
    fig_path = "figs/" + sub_dir
    if not path.isdir(fig_path):
        mkdir(fig_path)

    quantity = "|M|"
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.set_ylabel("$" + quantity + "/[\\mathrm{" + units[quantity] + "}]$")
    ax.set_xlabel("$T / [J]$")
    for i, N in enumerate(Ns):
        ax.plot(Ts, samples[quantity][i], styles[i], label="$N={}$".format(int(N)))
    Ts2 = np.linspace(1.5, Tc, 1000)
    ax.plot(Ts2, (1 - np.sinh(2/ Ts2)**(-4))**(1 / 8), "--", label="Analytical sol.")
    ax.legend()
    ax.grid(True)
    
    plt.savefig(fig_path + "analytic.png")
    plt.close(fig)

def time_dependence(sub_dir):
    pass