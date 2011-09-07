# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.title("Number of DOF as a function of physical time")
pylab.xlabel("Physical time")
pylab.ylabel("Number of DOF")
axis('equal')
data = numpy.loadtxt("conv_dof_est.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="hp-FEM")
legend()

# initialize new window
pylab.figure()

# plot CPU convergence graph
pylab.title("CPU time as a function of physical time")
pylab.xlabel("Physical time")
pylab.ylabel("CPU time")
axis('equal')
data = numpy.loadtxt("conv_cpu_est.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="hp-FEM")
legend()

# finalize
show()
