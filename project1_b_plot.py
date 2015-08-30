#
#     Project1: b)
#     The algorithm for solving the tridiagonal matrix
#     equation is implemented in project1_b.cpp.
#     Results are read from textfile and plotted here.
#
from math import *
import numpy as np
import matplotlib.pyplot as plt

def read_x_u_v(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    x = []; u = []; v = [];
    # Read lines except for the first one:
    lines = infile.readlines()[1:]
    for line in lines:
    	words = line.split()
        x.append(float(words[0]))
        u.append(float(words[1]))
        v.append(float(words[2]))
    infile.close()
    return x, u, v

# Fetching data by a call on read_x_u_v for three different n:
x1, u1, v1 = read_x_u_v('data_n10.txt')
x2, u2, v2 = read_x_u_v('data_n100.txt')
x3, u3, v3 = read_x_u_v('data_n1000.txt')


# Plotting commands:
plt.rcParams.update({'font.size': 16})
fig, ax = plt.subplots(1)
ax.plot(x1,u1,'k--',label='Analytical, $u(x)$')
ax.plot(x1,v1,'r-',label='Numerical, $v(x)$')
ax.set_xlabel('$x\in (0,1)$')
ax.set_ylabel('Solution of differential equation: $u(x)$')
ax.legend(loc='lower center',fancybox='True')
ax.set_title('Plot of numerical and analytical solution for $n = 10$')
ax.grid()
plt.show()

fig, ax = plt.subplots(1)
ax.plot(x2,u2,'k--',label='Analytical, $u(x)$')
ax.plot(x2,v2,'m-',label='Numerical, $v(x)$')
ax.set_xlabel('$x\in (0,1)$')
ax.set_ylabel('Solution of differential equation: $u(x)$')
ax.legend(loc='lower center',fancybox='True')
ax.set_title('Plot of numerical and analytical solution for $n = 100$')
ax.grid()
plt.show()

fig, ax = plt.subplots(1)
ax.plot(x3,u3,'k--',label='Analytical, $u(x)$')
ax.plot(x3,v3,'g-',label='Numerical, $v(x)$')
ax.set_xlabel('$x\in (0,1)$')
ax.set_ylabel('Solution of differential equation: $u(x)$')
ax.legend(loc='lower center',fancybox='True')
ax.set_title('Plot of numerical and analytical solution for $n = 1000$')
ax.grid()
plt.show()
