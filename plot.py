import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import sys
from textwrap import wrap
import matplotlib.patches as mpatches
from collections import defaultdict




def serial_draw1():
    f_a = open("serial_sum.txt")

    read_a = f_a.read()
    split_a = read_a.split('\n')

    # Get data out of the data file
    size_a = []
    time_a = []

    for item in split_a:
        if item:
            size_time_rate = item.split(' ')
            size_a.append(int(size_time_rate[0]))
            time_a.append(float(size_time_rate[1]))

    # Plot the graph
    title2 = ('log-log plot of running time vs the number of particles')
    plt.plot([0, 160000], [0, 24.44256], '-r', label='Ideal performance')
    plt.plot(size_a, time_a, '-g', label='Optimized serial')
    plt.plot(size_a, time_a, '^b')
    plt.legend()

    plt.xlabel('Particle Sizes')
    plt.ylabel('Running Time')

    plt.xscale('log')
    plt.yscale('log')

    plt.subplots_adjust(top=0.85)
    plt.title(wrap(title2, 60), y = 1.08)
    plt.savefig('serial_2.png')

    # Clean up
    plt.close()
    f_a.close()

#serial_draw1()

def omp_strong_time_draw():
    f_a = open("openmp_sum.txt")

    read_a = f_a.read()
    split_a = read_a.split('\n')

    # Get data out of the data file
    size_a = []
    time_a = []
    processor_a = []

    for item in split_a:
        if item:
            size_time_rate = item.split(' ')
            size_a.append(int(size_time_rate[0]))
            processor_a.append(int(size_time_rate[1]))
            time_a.append(float(size_time_rate[2]))

    time_2000 = []
    processor_2000 = []
    time_10000 = []
    processor_10000 = []
    for i in range(len(size_a)):
        if size_a[i] == 2000:
            time_2000.append(time_a[i])
            processor_2000.append(processor_a[i])
        elif size_a[i] == 10000:
            time_10000.append(time_a[i])
            processor_10000.append(processor_a[i])

    legends = []
    patch = mpatches.Patch(color='b', label='particle size: 2000')
    legends += [patch]
    patch = mpatches.Patch(color='g', label='particle size: 10000')
    legends += [patch]

    # Plot the graph
    # title2 = ('strong scaling of running time - communication pattern')
    title2 = ('communication time based on strong scaling experiment')

    plt.plot(processor_2000, time_2000, '-b')
    plt.plot(processor_2000, time_2000, '^b')

    plt.plot(processor_10000, time_10000, '-g')
    plt.plot(processor_10000, time_10000, '^g')

    plt.legend(handles=legends)

    plt.xlabel('Num of Processor')
    plt.ylabel('Running Time (second)')

    # plt.xscale('log')
    # plt.yscale('log')

    plt.subplots_adjust(top=0.85)
    plt.title(wrap(title2, 60), y = 1.08)
    plt.savefig('mpi_strong.png')

    # Clean up
    plt.close()
    f_a.close()

omp_strong_time_draw()

colours = {1 : 'r',
           2000 : 'b',
           4 : 'm',
           8 : 'k',
           10000: 'g'
           }


def omp_weak_time_draw():
    f_a = open("openmp_sum.txt")

    read_a = f_a.read()
    split_a = read_a.split('\n')

    # Get data out of the data file
    size_a = []
    time_a = []
    processor_a = []

    for item in split_a:
        if item:
            size_time_rate = item.split(' ')
            size_a.append(int(size_time_rate[0]))
            processor_a.append(int(size_time_rate[1]))
            time_a.append(float(size_time_rate[2]))

    times = defaultdict(list)
    sizes = defaultdict(list)

    for i in range(len(size_a)):
        if size_a[i]/processor_a[i] == 2000:
            sizes[2000].append(size_a[i])
            times[2000].append(time_a[i])
        elif size_a[i]/processor_a[i] == 10000:
            sizes[10000].append(size_a[i])
            times[10000].append(time_a[i])


    work_processor = [2000, 10000]

    legends = []
    for num in work_processor:
        patch = mpatches.Patch(color=colours[num], label='Work/Processor: %d' % num)
        legends += [patch]

    # Plot the graph
    #title2 = ('weak scaling of running time - communication pattern')
    title2 = ('communication time computed based on weak scaling experiment')

    plt.clf()
    for num in work_processor:
        plt.plot(sizes[num], times[num], colours[num])
        plt.plot(sizes[num], times[num], '^'+colours[num])

    plt.legend(handles=legends)

    plt.xlabel('Particle Sizes')
    plt.ylabel('Running Time (second)')

    plt.subplots_adjust(top=0.85)
    plt.title(wrap(title2, 60), y = 1.08)
    plt.savefig('mpi_weak.png')

    # Clean up
    plt.close()
    f_a.close()

omp_weak_time_draw()

def omp_strong_speed():
    # Get data out of the data file
    num_thread = [1, 2, 4, 8, 16]
    speedup_large = [0.79,    1.55,    2.79,    4.14,   4.73]
    speedup_small = [0.76,    1.43,    2.21,    2.43,   1.45]
    # Plot the graph
    title = ('Speedup VS num of threads')
    plt.plot(num_thread, speedup_large, '-g', label='particle size=10000')
    plt.plot(num_thread, speedup_large, '^g')
    plt.plot(num_thread, speedup_small, '-b', label='particle size=2000')
    plt.plot(num_thread, speedup_small, '^b')

    
    plt.legend()
    plt.xlabel('Num of threads')
    plt.ylabel('Speedup')

    plt.subplots_adjust(top=0.85)
    plt.title(wrap(title, 60), y = 1.08)
    plt.savefig('mpi_strong_speedup.png')

    # Clean up
    plt.close()

# omp_strong_speed()

def omp_efficiency():
    # Get data out of the data file
    num_thread = [1, 2, 4, 8, 16]

    strong_large = [0.79,    0.77,    0.70,    0.52,    0.30]
    weak_large = [0.79,    0.79,    0.77,    0.73,    0.62]

    strong_small = [0.76,    0.71,    0.55,    0.30,    0.09]
    weak_small = [0.76,    0.76,    0.72,    0.66,    0.52]


    # Plot the graph
    title = ('Strong scaling efficiency VS num of threads')

    plt.plot(num_thread, strong_large, '-g', label='particle size=10000')
    plt.plot(num_thread, strong_large, '^g')
    plt.plot(num_thread, strong_small, '-b', label='particle size=2000')
    plt.plot(num_thread, strong_small, '^b')
    plt.legend()
    plt.xlabel('Num of threads')
    plt.ylabel('Efficiency')

    plt.subplots_adjust(top=0.85)
    plt.title(wrap(title, 60), y = 1.08)
    plt.savefig('mpi_strong_efficiency.png')

    # Clean up
    plt.close()

    # Plot the graph
    title = ('Weak scaling efficiency VS num of threads')
    plt.plot(num_thread, weak_large, '-g', label='particle size=10000')
    plt.plot(num_thread, weak_large, '^g')
    plt.plot(num_thread, weak_small, '-b', label='particle size=2000') 
    plt.plot(num_thread, weak_small, '^b')
    plt.legend()
    plt.xlabel('Num of threads')
    plt.ylabel('Efficiency')

    plt.subplots_adjust(top=0.85)
    plt.title(wrap(title, 60), y = 1.08)
    plt.savefig('mpi_weak_efficiency.png')

    # Clean up
    plt.close()

# omp_efficiency()

# rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)
#
# threads, timings = np.loadtxt('strong_scaling.csv', delimiter=',', usecols=(0,1), unpack=True)
#
# serial_time = timings[0];
# timings = np.divide(serial_time, timings)
#
# plt.plot(threads, timings,'k')
# plt.xlim([1,26])
# plt.xlabel("Number of OMP threads")
# plt.ylabel("Speedup over the serial implementation")
# plt.savefig('strong_scaling.png', dpi=300)
# plt.clf()
# plt.cla()
#
# threads, timings = np.loadtxt('weak_scaling.csv', delimiter=',', usecols=(0,1), unpack=True)
# serial_time = timings[0];
# timings = np.divide(serial_time, timings)
#
# plt.plot(threads, timings,'k')
# plt.xlim([1,20])
# plt.xlabel("Number of OMP threads")
# plt.ylabel("Efficiency")
# plt.savefig('weak_scaling.png', dpi=300)
# plt.clf()
# plt.cla()