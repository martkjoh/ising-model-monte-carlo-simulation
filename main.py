from utilities import get_samples, write_samples
from progress.bar import Bar
from observables import *  


equib = 100
n = 100
m = 20
Ns = [8, 16]
Ts = np.linspace(1.5, 1.2 * Tc, m)

observables = {
    "E": energy_density,
    "E2": energy_density_squared,
    "abs_M": absolute_magnetization,
    "M2": magnetization_squared,
}

def sample_observables():
    """ creates a file with samples from a MC-walk of the observables """
    bar = Bar("sampling", max=len(Ns)*m)
    samples = {key : np.empty(m) for key in observables}
    
    for _, N in enumerate(Ns):
        for i, T in enumerate(Ts):
            data = get_samples(N, T, n, equib, observables)
            for j, key in enumerate(samples):
                samples[key][i] = data[j]
        bar.next()
    
    bar.finish()
    print("writing data")
    write_samples(samples, Ns, Ts)