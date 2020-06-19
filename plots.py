import numpy as np
from matplotlib import pyplot as plt
from os import path, mkdir
from progress.bar import Bar

from utilities import get_s, MC_sweep, get_samples, read_samples, read_times, read
from physical_quantities import observables, func_of_obs, Tc, units, latex, reudced_temp, tension


font = {
    'family': 'serif',
    'weight': 'normal',
    'size': 30}
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rc("lines", linewidth=1, markersize=3)
styles = ["--x", "--o", "--+", "--d", "--1", "--*"]


def equilibration():
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


def vals(sub_dir):
    Ns, Ts, samples = read_samples(list(observables), sub_dir=sub_dir)
    fig_path = "figs/" + sub_dir
    if not path.isdir(fig_path):
        mkdir(fig_path)

    for quantity in list(observables):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_ylabel("$" + latex[quantity] + "/[\\mathrm{" + units[quantity] + "}]$")
        ax.set_xlabel("$T / [J]$")

        for i, N in enumerate(Ns):
            ax.plot(Ts, samples[quantity][i] / 2, styles[i%len(styles)], label="$N={}$".format(int(N)))
        ax.plot([Tc, Tc], ax.get_ylim(), "k--", label="$T_c$")
        ax.grid(True)
        ax.legend()        
        
        plt.tight_layout()
        plt.savefig(fig_path + quantity + ".png", dpi=300)
        plt.close(fig)

def funcs(sub_dir=""):
    Ns, Ts, samples = read_samples(list(observables), sub_dir=sub_dir)
    fig_path = "figs/" + sub_dir
    if not path.isdir(fig_path):
        mkdir(fig_path)

    for quantity in list(func_of_obs):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_ylabel("$" + latex[quantity] + "/[\\mathrm{" + units[quantity] + "}]$")
        ax.set_xlabel("$T / [J]$")

        y = func_of_obs[quantity](samples, Ts, Ns)
        for i, N in enumerate(Ns):
            ax.plot(Ts, y[i] / 2, styles[i%len(styles)], label="$N={}$".format(int(N)))
        ax.plot([Tc, Tc], ax.get_ylim(), "k--", label="$T_c$")
        ax.grid(True)
        ax.legend()

        plt.tight_layout()
        plt.savefig(fig_path + quantity + ".png", dpi=300)
        plt.close(fig)
        

def susc(sub_dir):
    Ns, Ts, samples = read_samples(list(observables), sub_dir=sub_dir)
    fig_path = "figs/" + sub_dir
    if not path.isdir(fig_path):
        mkdir(fig_path)

    quantity = "absolute_magnetization"
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.set_ylabel("$" + quantity + "/[\\mathrm{" + units[quantity] + "}]$")
    ax.set_xlabel("$T / [J]$")
    for i, N in enumerate(Ns):
        ax.plot(Ts, samples[quantity][i], styles[i%len(styles)], label="$N={}$".format(int(N)))
    Ts2 = np.linspace(1.5, Tc, 1000)
    ax.plot(Ts2, (1 - np.sinh(2/ Ts2)**(-4))**(1 / 8), "--", label="Analytical sol.")
    ax.legend()
    ax.grid(True)
    
    plt.savefig(fig_path + "analytic.png")
    plt.close(fig)


def time_dependence(sub_dir):
    times, times_MJ = read_times(sub_dir)
    Ns, _, _ = read_samples(list(observables), sub_dir=sub_dir)
    fig_path = "figs/" + sub_dir

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(Ns, times, "x", label="$t$")
    ax.plot(Ns, times_MJ, "+", label="$t$")
    plt.savefig(fig_path + "times.png")
    plt.close(fig)


def Mon_Jasnow(sub_dir):
    Ns = read(sub_dir, "sizes")
    Ts = read(sub_dir, "temps")
    name = "tau_MJ"
    tau = read(sub_dir, name)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.set_ylabel("$\\tau/[J]$")
    ax.set_xlabel("$T/[J]$")
    for i, N in enumerate(Ns):
        ax.plot(Ts, tau[i], "--.", label="$N={}$".format(int(N)))
    ax.plot([Tc, Tc], ax.get_ylim(), "k--", label="$T_c$")
    ax.legend()
    ax.grid(True)
    
    plt.tight_layout()
    plt.savefig("figs/" + sub_dir + name + ".png", dpi=300)
    plt.close(fig)



def more(sub_dir):
    """ What the fuuuuck??!?! """

    Ns, Ts, samples = read_samples(list(observables), sub_dir=sub_dir)
    name = "tau_MJ"
    tau = tension(samples, Ts, Ns)
    # tau = read(sub_dir, name)
    fig_path = "figs/" + sub_dir
    if not path.isdir(fig_path):
        mkdir(fig_path)
    
    fig, ax = plt.subplots()
    m = 38
    ax.set_title("$T = {}$".format(Ts[m]))
    ax.loglog(tau[:, m], Ns)
    plt.show()

    fig, ax = plt.subplots()
    for i, t in enumerate(reudced_temp(Ts[30:40:2])):
        ax.plot(tau[:, i] / t, 1 / (Ns * t), "--+", label="$t = {}$".format(t))
        ax.legend()
    plt.show()

    fig, ax = plt.subplots(2)
    print(np.shape(tau))
    for i, N in enumerate(Ns):
        ax[0].plot(Ts, tau[i], styles[i], label="$N={}$".format(int(N)))
    for i, T in enumerate(Ts):
        ax[1].plot(Ns, Ns*tau[:, i], "--x", label="$T={}$".format(int(T)))
    plt.show()


if __name__ == "__main__":
    more("quantities/")
