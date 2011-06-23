Screen (Maxwell's Equations)
----------------------------

**Git reference:** Benchmark `screen <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/benchmarks-general/screen>`_.

This example solves time-harmonic Maxwell's equations. It describes an electromagnetic wave that 
hits a thin screen under the angle of 45 degrees, causing a singularity at the tip of the screen.
The strength of the singularity makes this example rather difficult. 

Model problem
~~~~~~~~~~~~~

Equation solved: Time-harmonic Maxwell's equations

.. math::
    :label: screen

    \frac{1}{\mu_r} \nabla \times \nabla \times E - \kappa^2 \epsilon_r E = \Phi.

Domain of interest is the square $(-1,1)^2$ missing the edge that connects the center with 
the midpoint of the left side. It is filled with air:

.. image:: benchmark-screen/domain.png
   :align: center
   :width: 490
   :alt: Computational domain.

Boundary conditions
~~~~~~~~~~~~~~~~~~~

Tangential component of solution taken from known exact solution (essential BC). 

Exact solution 
~~~~~~~~~~~~~~

This is rather complicated in this case - see the file 
`definitions.cpp <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/benchmarks-general/screen/definitions.cpp>`_.

Sample solution
~~~~~~~~~~~~~~~

Real part of $E_1$:

.. image:: benchmark-screen/sol1.png
   :align: center
   :width: 510
   :alt: Solution.

Real part of $E_2$:

.. image:: benchmark-screen/sol2.png
   :align: center
   :width: 510
   :alt: Solution.

Imaginary part of $E_1$:

.. image:: benchmark-screen/sol3.png
   :align: center
   :width: 510
   :alt: Solution.

Imaginary part of $E_2$:

.. image:: benchmark-screen/sol4.png
   :align: center
   :width: 510
   :alt: Solution.

Convergence comparisons
~~~~~~~~~~~~~~~~~~~~~~~

Final mesh (h-FEM with linear elements):

.. image:: benchmark-screen/mesh-h1.png
   :align: center
   :width: 460
   :alt: Final mesh (h-FEM with linear elements).

Note that the polynomial order indicated corresponds to the tangential components 
of approximation on element interfaces, not to polynomial degrees inside the elements
(those are one higher).

Final mesh (h-FEM with quadratic elements):

.. image:: benchmark-screen/mesh-h2.png
   :align: center
   :width: 460
   :alt: Final mesh (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: benchmark-screen/mesh-hp.png
   :align: center
   :width: 460
   :alt: Final mesh (hp-FEM).

DOF convergence graphs:

.. image:: benchmark-screen/conv_dof.png
   :align: center
   :width: 600
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: benchmark-screen/conv_cpu.png
   :align: center
   :width: 600
   :alt: CPU convergence graph.
   
   
