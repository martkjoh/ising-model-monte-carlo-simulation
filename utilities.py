import numpy as np
from numpy import exp, sqrt, log as ln
from numpy.random import choice, randint, rand
from physical_quantities import ms, energy_MJ
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


def get_delta_H_MJ(s, N):
    """ retrns lattice with change in energy from flipping each spin"""

    sum_neigh = np.zeros_like(s)
    s_pp = np.concatenate([np.ones((1, N), dtype=int), s, np.ones((1, N), dtype=int)])
    for j in range(2):
        for n in (-1, 1):
            sum_neigh += np.roll(s_pp, n, axis=j)[1: -1]

    return 2*s*sum_neigh


def MC_sweep(s, N, T, dH=get_delta_H):
    """ MC-hastings sweep, with prob. c of trying to flipping each spin each loop """

    c = 0.5
    for _ in range(int(1 / c)):
        to_flip = rand(N, N) < c
        delta_H = dH(s, N)
        pos = delta_H > 0
        to_flip2 = exp(-delta_H[pos] / T) > rand(np.sum(pos))
        to_flip[pos] = np.logical_and(to_flip[pos], to_flip2)
        s[to_flip] *= -1


def get_sample_from_state(s, N, T, observables):
    """ samples a set of observables from a configuration of s """

    sample = np.empty(len(observables))
    for i, key in enumerate(observables):
        sample[i] = observables[key](s, N, T)
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
        sample_avg += get_sample_from_state(s, N, T, observables)
    return sample_avg / n


def Mon_Jasnow(N, T, n, equib):
    """
    Runs the Mon-Jasnow algorithm on an N*N matrix at temperature T, 
    first equilibrating it for equib MC-sweeps, then sampling after each of n steps.
    """

    fraction = 0

    s = get_s(N)
    for _ in range(equib):
        MC_sweep(s, N, T, dH=get_delta_H_MJ)
    for _ in range(n):
        MC_sweep(s, N, T, dH=get_delta_H_MJ)
        fraction += exp(-2 / T * ms(s, N))
    return - T / N * ln(fraction / n)


def write(data, sub_dir, name):
    data_path = "data/" + sub_dir

    if not path.isdir(data_path):
        mkdir(data_path)
    np.savetxt(data_path + name + ".csv", data)


def read(sub_dir, name):
    data_path = "data/" + sub_dir
    if not path.isdir(data_path):
        raise  Exception("Direcotry does not exist")
    return np.loadtxt(data_path + name + ".csv")


def write_samples(samples, Ns, Ts, times, sub_dir):
    """ Writes a dict of samples at different Ns and Ts to .csv files """

    write(Ns, sub_dir, "sizes")
    write(Ts, sub_dir, "temps")
    write(times, sub_dir, "times")
    for key in samples:
        write(samples[key], sub_dir, key)


def read_samples(data, sub_dir):
    """ Reads a dict of samples at different Ns and Ts from .csv files """

    Ns = read(sub_dir, "sizes")
    Ts = read(sub_dir, "temps")
    samples = {name : 0 for name in data}
    for name in data:
        samples[name] = read(sub_dir, name)
        
    return Ns, Ts, samples


def read_times(sub_dir):
    data_path = "data/" + sub_dir

    if not path.isdir(data_path):
        raise Exception("Direcotry does not exist")
    
    return np.loadtxt(data_path + "times.csv"), np.loadtxt(data_path + "times_MJ.csv")