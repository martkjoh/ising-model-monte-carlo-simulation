from utilities import get_samples, write_samples, write, Mon_Jasnow
import plots
from progress.bar import Bar
import numpy as np
from physical_quantities import observables, Tc
from time import time

equib = 10_000
n = 1_000
sub_dir = "quantities/"

# number of different tempratures to simulate
temps = 20
Ts = np.linspace(1.5, 1.2*Tc, temps)
# The different sizes of the grid to simutale
Ns = [8, 16, 32, 64]
sizes = len(Ns)


def sample_observables(sub_dir):
    """ creates a file with samples of the observables from a MC-walk """

    bar = Bar("sampling", max=sizes * temps)
    samples_dict = {key: np.empty((sizes, temps)) for key in observables}
    times = np.zeros(sizes)
    
    for i, N in enumerate(Ns):
        for j, T in enumerate(Ts):
            t = time()

            samples = get_samples(N, T, n, equib, observables)
            for k, key in enumerate(samples_dict):
                samples_dict[key][i, j] = samples[k]
            
            bar.next()
            times[i] += time() - t

    
    bar.finish()
    print("writing data")
    times /= temps
    write_samples(samples_dict, Ns, Ts, times, sub_dir=sub_dir)


def get_tension(sub_dir):
    tau = np.empty((sizes, temps))
    bar = Bar(max=temps*sizes)
    for i, N, in enumerate(Ns):
        for j, T in enumerate(Ts):
            tau[i, j] =  Mon_Jasnow(N, T, n, equib)
            bar.next()
    bar.finish()

    write(tau, sub_dir, "tau")


def full_suite(sub_dir="test/", gen_data=False):
    """ Runs the full suit of functinos: generating data and making plots """

    if gen_data:
        # sample_observables(sub_dir)
        get_tension(sub_dir)
    
    plots.vals(sub_dir)
    plots.funcs(sub_dir)
    plots.susc(sub_dir)
    plots.time_dependence(sub_dir)
    plots.Mon_Jasnow(sub_dir)


if __name__ == "__main__":
    # plot_equilibration()

    full_suite(gen_data=True)
