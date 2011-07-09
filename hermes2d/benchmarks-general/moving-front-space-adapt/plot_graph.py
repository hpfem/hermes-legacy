# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.title("DOF history")
pylab.xlabel("Physical time (s)")
pylab.ylabel("Number of degrees of freedom")
data = numpy.loadtxt("dof_history.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, "-s", label="number of DOF")

legend()

# finalize
show()
