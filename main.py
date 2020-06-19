from utilities import get_samples, write_samples, write, Mon_Jasnow
import plots
from progress.bar import Bar
import numpy as np
from physical_quantities import observables, Tc
from time import time

equib = 100_000
n = 100_000
sub_dir = "quantities/"

# number of different tempratures to simulate
temps = 50
Ts = np.linspace(0.5, 1.2*Tc, temps)
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
    write_samples(samples_dict, Ns, Ts, times, sub_dir=sub_dir)


def run_mon_jas(sub_dir):
    tau = np.empty((sizes, temps))
    bar = Bar("sampling", max=temps*sizes)
    for i, N, in enumerate(Ns):
        for j, T in enumerate(Ts):
            tau[i, j] =  Mon_Jasnow(N, T, n, equib)
            bar.next()
    bar.finish()

    write(tau, sub_dir, "tau_MJ")


def full_suite(sub_dir="test/", gen_data=False):
    """ Runs the full suit of functinos: generating data and making plots """

    if gen_data:
        sample_observables(sub_dir)
        run_mon_jas(sub_dir)
    
    plots.vals(sub_dir)
    plots.funcs(sub_dir)
    plots.susc(sub_dir)
    plots.time_dependence(sub_dir)
    plots.Mon_Jasnow(sub_dir)


if __name__ == "__main__":
    full_suite(sub_dir = sub_dirgen_data=True)
