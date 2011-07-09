# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.title("Error history")
pylab.xlabel("Physical time")
pylab.ylabel("Error")
data = numpy.loadtxt("err_est_history.dat")
x = data[:, 0]
y = data[:, 1]
semilogy(x, y, "-", label="error estimate")
data = numpy.loadtxt("err_exact_history.dat")
x = data[:, 0]
y = data[:, 1]
semilogy(x, y, "-", label="exact error")

legend()

# initialize new window
#pylab.figure()



# finalize
show()
