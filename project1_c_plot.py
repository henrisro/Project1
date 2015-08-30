#
#     Project1: c)
#     The algorithm for solving the tridiagonal matrix
#     equation is implemented in project1_c.cpp.
#     Log of the absolute error is read from file and plotted:
#
from math import *
import numpy as np
import matplotlib.pyplot as plt

def read_logh_logeps(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    log_h = []; log_eps = [];
    # Read lines except for the first one:
    lines = infile.readlines()[1:]
    for line in lines:
    	words = line.split()
        log_h.append(float(words[1]))
        log_eps.append(float(words[2]))
    infile.close()
    return log_h, log_eps

# Fetching data and plotting them:
log_h, log_eps = read_logh_logeps('log_error.txt')

# Plotting commands:
plt.rcParams.update({'font.size': 16})
fig, ax = plt.subplots(1)
ax.plot(log_h,log_eps,'bo-')
ax.set_xlabel('$\log_{10}(h)$')
ax.set_ylabel('max$(\log_{10}(\mid (u_i-v_i)/u_i \mid ))$')
ax.set_title('Plot of maximal relative error as function of step size')
ax.grid()
plt.show()