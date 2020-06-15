import numpy as np
from numpy import exp, sqrt, log as ln
from numpy.random import choice, randint, rand


def get_s(N):
    """ Return lattice of spins, randomly oriented """
    return choice((-1, 1), (N, N))

def get_delta_H(s, N):
    """ retrns lattice with change in energy from flipping each spin"""
    sum_neigh = np.zeros_like(s)
    for j in range(2):
        for n in (-1, 1):
            sum_neigh += np.roll(s, n, axis=j)

    return 2*s*sum_neigh

def MC_sweep(s, N, T):
    """ MC-hastings sweep, with prob. c of trying to flipping each spin each loop """
    c = 0.5
    for _ in range(int(1 / c)):
        to_flip = rand(N, N) < c
        delta_H = get_delta_H(s, N)
        pos = delta_H > 0
        to_flip2 = exp(-delta_H[pos] / T) > rand(np.sum(pos))
        to_flip[pos] = np.logical_and(to_flip[pos], to_flip2)
        s[to_flip] *= -1


def sample(s, N, observables):
    """ samples a set of observables from a configuration of s """
    obs = np.empty(len(observables))
    for i, key in enumerate(observables):
        obs[i] = observables[key](s, N)
    return obs


def get_samples(N, T, n, equib, observables):
    """ 
    Creates a N*N matrix, equilibrates them for equib steps, then does
    n MC-sweeps and at temprature T samples all the observables, 
    returning an array of the samples
    """
    obs = np.zeros(len(observables))

    s = get_s(N)
    for _ in range(equib):
        MC_sweep(s, N, T)
    for _ in range(n):
        MC_sweep(s, N, T)
        obs += sample(s, N, observables)
    return obs / n

def write_samples(samples, Ns, Ts):
    """ Writes a dict of samples at different Ns and Ts to .csv files """
    for i, key in samples:
        np.savetxt("data" + key + ".csv", [Ns[i], Ts[i], samples[key]])