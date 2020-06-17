from utilities import get_samples, write_samples
from progress.bar import Bar
import numpy as np
from physical_quantities import observables, Tc
from plots import *
from time import time

equib = 100_000
n = 100_000
sub_dir = "quantities/"

# number of different tempratures to simulate
temps = 50
Ts = np.linspace(1.5, 1.5*Tc, temps)
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


def full_suite(sub_dir="test/", gen_data=False):
    """ Runs the full suit of functinos: generating data and making plots """

    if gen_data:
        sample_observables(sub_dir)
    
    plot_observables(sub_dir)
    plot_funcs(sub_dir)
    mag_plot(sub_dir)
    plot_time_dependence(sub_dir)


if __name__ == "__main__":
    # plot_equilibration()

    full_suite(sub_dir=sub_dir)
