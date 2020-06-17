import numpy as np
from numpy import exp, sqrt, log as ln
from numpy.random import choice, randint, rand
from os import path, mkdir

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


def get_sample_from_state(s, N, observables):
    """ samples a set of observables from a configuration of s """

    sample = np.empty(len(observables))
    for i, key in enumerate(observables):
        sample[i] = observables[key](s, N)
    return sample


def get_samples(N, T, n, equib, observables):
    """ 
    Creates a N*N matrix, equilibrates them for equib steps, then does
    n MC-sweeps and at temprature T samples all the observables, 
    returning an array of the samples
    """

    sample_avg = np.zeros(len(observables))
    s = get_s(N)
    for _ in range(equib):
        MC_sweep(s, N, T)
    for _ in range(n):
        MC_sweep(s, N, T)
        sample_avg += get_sample_from_state(s, N, observables)
    return sample_avg / n


def write_samples(samples, Ns, Ts, times, sub_dir):
    """ Writes a dict of samples at different Ns and Ts to .csv files """

    data_path = "data/" + sub_dir

    if not path.isdir(data_path):
        mkdir(data_path)

    np.savetxt(data_path + "sizes.csv", Ns)
    np.savetxt(data_path + "temps.csv", Ts)
    np.savetxt(data_path + "times.csv", times)
    for key in samples:
        np.savetxt(data_path + key + ".cvs", samples[key])


def read_samples(data, sub_dir):
    """ Reads a dict of samples at different Ns and Ts from .csv files """

    data_path = "data/" + sub_dir

    if not path.isdir(data_path):
        raise  Exception("Direcotry does not exist")

    samples = {name : 0 for name in data}

    Ns = np.loadtxt(data_path + "sizes.csv")
    Ts = np.loadtxt(data_path + "temps.csv")
    for name in data:
        samples[name] = np.loadtxt(data_path + name + ".cvs")
        
    return Ns, Ts, samples

def read_times(sub_dir):
    data_path = "data/" + sub_dir

    if not path.isdir(data_path):
        raise Exception("Direcotry does not exist")
    
    return np.loadtxt(data_path + "times.csv")