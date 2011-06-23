Interior Layer (Elliptic)
-------------------------

**Git reference:** Benchmark `layer-interior <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/benchmarks-general/layer-interior>`_.

This example has a smooth solution that exhibits a steep interior layer.

Model problem
~~~~~~~~~~~~~

Equation solved: Poisson equation 

.. math::
    :label: layer-interior

       -\Delta u - f = 0.

Domain of interest: Unit square $(0, 1)^2$.

Exact solution
~~~~~~~~~~~~~~

.. math::
    :label: layer-interior-exact

    u(x, y) = \mbox{atan}\left(S \sqrt{(x-1.25)^2 + (y+0.25)^2} - \pi/3\right).

where $S$ is a parameter (slope of the layer). With larger $S$, this problem 
becomes difficult for adaptive algorithms, and at the same time the advantage of 
adaptive $hp$-FEM over adaptive low-order FEM becomes more significant. We will 
use $S = 60$ in the following.

Right-hand side
~~~~~~~~~~~~~~~

Obtained by inserting the exact solution into the equation:

.. math::
    :label: layer-interior-rhs
 
    f(x, y) = \frac{27}{2} (2y + 0.5)^2 (\pi - 3t) \frac{S^3}{u^2 t_2} +
    \frac{27}{2} (2x - 2.5)^2 (\pi - 3t) \frac{S^3}{u^2 t_2}
    - \frac{9}{4} (2y + 0.5)^2 \frac{S}{u t^3} -
    \frac{9}{4} (2x - 2.5)^2 \frac{S}{u t^3} +
    18 \frac{S}{ut}.

Boundary conditions
~~~~~~~~~~~~~~~~~~~

Nonconstant Dirichlet, matching the exact solution.

Sample solution
~~~~~~~~~~~~~~~

.. image:: benchmark-layer-interior/sol_3d_view.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution.

Convergence comparison
~~~~~~~~~~~~~~~~~~~~~~

Final mesh (h-FEM with linear elements):

.. image:: benchmark-layer-interior/mesh-h1.png
   :align: center
   :width: 500
   :alt: Final mesh (h-FEM with linear elements).

Final mesh (h-FEM with quadratic elements):

.. image:: benchmark-layer-interior/mesh-h2.png
   :align: center
   :width: 500
   :alt: Final mesh (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: benchmark-layer-interior/mesh-hp.png
   :align: center
   :width: 500
   :alt: Final mesh (hp-FEM).

DOF convergence graphs:

.. image:: benchmark-layer-interior/conv_dof.png
   :align: center
   :width: 600
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: benchmark-layer-interior/conv_cpu.png
   :align: center
   :width: 600
   :alt: CPU convergence graph.
