from utilities import get_samples, write_samples
from progress.bar import Bar
from physical_quantities import observables
# from physical_quantities import observables
from plots import *
from time import time

equib = 100_000
n = 100_000
sub_dir = "quantities/"

# number of different tempratures to simulate
temps = 50
Ts = np.linspace(1.5, 1.5 * Tc, temps)
# The different sizes of the grid to simutale
Ns = [8, 16, 32, 64]
sizes = len(Ns)


def sample_observables(sub_dir=""):
    """ creates a file with samples of the observables from a MC-walk """

    bar = Bar("sampling", max=sizes * temps)
    samples_dict = {key: np.empty((sizes, temps)) for key in observables}
    times = np.zeros(sizes)
    
    for i, N in enumerate(Ns):
        for j, T in enumerate(Ts):
            t = time()
            samples = get_samples(N, T, n, equib, observables)
            times[i] += time() - t
            for k, key in enumerate(samples_dict):
                samples_dict[key][i, j] = samples[k]

            bar.next()
    
    bar.finish()
    print("writing data")
    times /= temps
    write_samples(samples_dict, Ns, Ts, times, sub_dir=sub_dir)


if __name__ == "__main__":
    # plot_equilibration()
    sample_observables(sub_dir=sub_dir)
    plot_observables(sub_dir=sub_dir)