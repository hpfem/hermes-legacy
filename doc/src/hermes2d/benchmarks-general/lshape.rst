L-Shape (Elliptic)
------------------

**Git reference:** Benchmark `lshape <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/benchmarks-general/lshape>`_.

This is a standard adaptivity benchmark whose exact solution is smooth but
contains singular gradient in a re-entrant corner. 

Model problem
~~~~~~~~~~~~~

Equation solved: Laplace equation 

.. math::
    :label: lshape

       -\Delta u = 0.

Domain of interest:

.. image:: benchmark-lshape/domain.png
   :align: center
   :width: 470
   :alt: Computational domain.

Exact solution
~~~~~~~~~~~~~~

.. math::
    :label: lshape-exact

    u(x, y) = r^{2/3}\sin(2a/3 + \pi/3)

where $r(x,y) = \sqrt{x^2 + y^2}$ and $a(x,y) = \mbox{atan}(x/y)$. 

Boundary conditions
~~~~~~~~~~~~~~~~~~~

Nonconstant Dirichlet, matching the exact solution.

Sample solution
~~~~~~~~~~~~~~~

.. image:: benchmark-lshape/sol_3d_view.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution.

Convergence comparisons
~~~~~~~~~~~~~~~~~~~~~~~

Final mesh (h-FEM with linear elements):

.. image:: benchmark-lshape/mesh-h1.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh (h-FEM with linear elements).

Final mesh (h-FEM with quadratic elements):

.. image:: benchmark-lshape/mesh-h2.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: benchmark-lshape/mesh-hp.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh (hp-FEM).

DOF convergence graphs:

.. image:: benchmark-lshape/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: benchmark-lshape/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.
