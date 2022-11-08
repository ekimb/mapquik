# get plots of # reads and % of correct reads per hit count

# example usage:
# first use get_q0_hitstats.sh [k] for all k,
# make sure mapquik-default-[k].paf and mapquik-default-[k]-Q0-incorrect.stats are in directory,
# python get_num_hits_per_read.py 5 6 7 8

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def read_incorrect_reads(path, read_incorrect):    
    with open(path) as fp:
        for line in fp:
            if '\t' in line:
                (read_id, hit_count) = line.split('\t')
                read_incorrect[read_id] = int(hit_count)

def read_all(path, dict_incorr):
    hit_all = {}
    with open(path) as fm:
        for line in fm:
            spl = line.split('\t')
            if len(spl) < 12: continue
            read_id = spl[0]
            hit_count = spl[9]
            if read_id in dict_incorr.keys():  # read is wrongly mapped
                if int(hit_count) in hit_all.keys():
                    (corr, incorr) = hit_all[int(hit_count)]
                    hit_all[int(hit_count)] = (corr, incorr + 1)
                else:
                    hit_all[int(hit_count)] = (1, 0)
            else:
                if int(hit_count) in hit_all.keys():
                    (corr, incorr) = hit_all[int(hit_count)]
                    hit_all[int(hit_count)] = (corr + 1, incorr)
                else:
                    hit_all[int(hit_count)] = (1, 0)
    return hit_all

def get_counts(num_correct_per_hit, num_incorrect_per_hit):
    total_hit_counts_set = set(list(num_correct_per_hit.keys()) + list(num_incorrect_per_hit.keys()))
    total_hit_counts = sorted(list(total_hit_counts_set))
    corr_incorr_per_hit_count = {}
    for hit_count in total_hit_counts:
        corr_count = 0
        incorr_count = 0
        if int(hit_count) in num_correct_per_hit.keys():
            corr_count = int(num_correct_per_hit[int(hit_count)])
        if int(hit_count) in num_incorrect_per_hit.keys():
            incorr_count = int(num_incorrect_per_hit[int(hit_count)])
        print("Hit count: ", hit_count, ": # correct: ", corr_count, " # incorrect: ", incorr_count)
        corr_incorr_per_hit_count[hit_count] = (corr_count, incorr_count)
    return corr_incorr_per_hit_count

# plot results

def plotter_percent(hit_all):
    x_ax = sorted([int(count) for count in hit_all.keys()])
    y_ax = []
    for count in x_ax:
        t = hit_all[count]
        y_ax.append((float(t[0]) / float(t[0] + t[1])) * 100.0)
    first_non100 = 0
    for i in range(len(y_ax) - 1, -1, -1):
        if y_ax[i] != 100.0:
            first_non100 = x_ax[i]
            break
            

    return x_ax, y_ax, first_non100

def plotter_reads(hit_all):
    x_ax = sorted([int(count) for count in hit_all.keys()])
    y_ax = []
    for count in x_ax:
        t = hit_all[count]
        y_ax.append(t[0] + t[1])
    return x_ax, y_ax

   
def plotter_up_to_100(hit_all, first_non100):
    x_ax = sorted([int(count) for count in hit_all.keys() if int(count) <= first_non100])
    y_ax = []
    for count in x_ax:
        t = hit_all[count]
        y_ax.append((float(t[0]) / float(t[0] + t[1])) * 100.0)
    return x_ax, y_ax


def obtain_stats_for_k(k):
    dict_incorr = {}
    incorr = "mapquik-default-" + k + "-Q0-incorrect.stats"
    allmaps = "mapquik-default-" + k + ".paf"
    read_incorrect_reads(incorr, dict_incorr)
    hit_all = read_all(allmaps, dict_incorr)
    x_ax_r, y_ax_r = plotter_reads(hit_all)
    x_ax_p, y_ax_p, first_non100 = plotter_percent(hit_all)
    x_ax_n, y_ax_n = plotter_up_to_100(hit_all, first_non100)
    return x_ax_r, y_ax_r, x_ax_p, y_ax_p, first_non100, x_ax_n, y_ax_n

def obtain_stats_for_all(axr, axn, k_list):
    axr.set_xlabel('Maximal pseudo-chain score per read')  
    axr.set_ylabel('# of reads')
    axn.set_xlabel('Maximal pseudo-chain score per read')  # Add an x-label to the axes.
    axn.set_ylabel('% of correctly-mapped reads')
    for k in k_list:
        x_ax_r, y_ax_r, x_ax_p, y_ax_p, first_non100, x_ax_n, y_ax_n = obtain_stats_for_k(k)
        axr.plot(x_ax_r, y_ax_r, label = "k = " + k)
        axn.plot(x_ax_n, y_ax_n, label = "k = " + k)
    
    axr.legend()
    axn.legend()
        


    






if len(sys.argv) == 1: print("Enter at least one k value."); exit;
k_list = sys.argv[1:]
#fig, (axr, axp) = plt.subplots(1, 2, figsize=(5, 2.7))
fig, (axr, axn) = plt.subplots(1, 2, figsize=(5, 2.7))
obtain_stats_for_all(axr, axn, k_list)
plt.show()

