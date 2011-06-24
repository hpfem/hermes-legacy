Saphir
------

**Git reference:** Example `saphir <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/examples/neutronics/saphir>`_.

Problem description
~~~~~~~~~~~~~~~~~~~

This is a standard nuclear engineering benchmark (IAEA number EIR-2) describing 
an external-force-driven configuration without fissile materials present, using one-group 
neutron diffusion approximation

.. math::
    :label: saphir

       -\nabla \cdot (D(x,y) \nabla \Phi) + \Sigma_a(x,y) \Phi - Q_{ext}(x,y) = 0.

The domain of interest is a 96 x 86 cm rectangle consisting of five regions:

.. image:: example-saphir/saphir.png
   :align: center
   :width: 400
   :alt: Schematic picture for the saphir example.

The unknown is the neutron flux $\Phi(x, y)$. The values of the diffusion coefficient 
$D(x, y)$, absorption cross-section $\Sigma_a(x, y)$ and the source term $Q_{ext}(x,y)$
are constant in the subdomains. The source $Q_{ext} = 1$ in areas 1 and 3 and zero 
elsewhere. Boundary conditions for the flux $\Phi$ are zero everywhere. 

This example uses multiple weak forms that are associated with different material 
markers.

Sample results
~~~~~~~~~~~~~~

Solution:

.. image:: example-saphir/saphir-sol.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution to the saphir example.

Final mesh (h-FEM with linear elements):

.. image:: example-saphir/saphir-mesh-h1.png
   :align: center
   :width: 440
   :alt: Final finite element mesh for the saphir example (h-FEM with linear elements).

Final mesh (h-FEM with quadratic elements):

.. image:: example-saphir/saphir-mesh-h2.png
   :align: center
   :width: 440
   :alt: Final finite element mesh for the saphir example (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: example-saphir/saphir-mesh-hp.png
   :align: center
   :width: 440
   :alt: Final finite element mesh for the saphir example (hp-FEM).

DOF convergence graphs:

.. image:: example-saphir/conv_dof.png
   :align: center
   :width: 600
   :alt: DOF convergence graph for example saphir.

CPU time convergence graphs:

.. image:: example-saphir/conv_cpu.png
   :align: center
   :width: 600
   :alt: CPU convergence graph for example saphir.

