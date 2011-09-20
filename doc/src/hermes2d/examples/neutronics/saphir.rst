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

.. figure:: example-saphir/saphir.png
   :align: center
   :scale: 35%
   :figclass: align-center
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

.. figure:: example-saphir/saphir-sol.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Solution to the saphir example.

Final mesh (h-FEM with linear elements):

.. figure:: example-saphir/saphir-mesh-h1.png
   :align: center
   :scale: 40% 
   :figclass: align-center
   :alt: Final finite element mesh for the saphir example (h-FEM with linear elements).

Final mesh (h-FEM with quadratic elements):

.. figure:: example-saphir/saphir-mesh-h2.png
   :align: center
   :scale: 40% 
   :figclass: align-center
   :alt: Final finite element mesh for the saphir example (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. figure:: example-saphir/saphir-mesh-hp.png
   :align: center
   :scale: 40% 
   :figclass: align-center
   :alt: Final finite element mesh for the saphir example (hp-FEM).

DOF convergence graphs:

.. figure:: example-saphir/conv_dof.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: DOF convergence graph for example saphir.

CPU time convergence graphs:

.. figure:: example-saphir/conv_cpu.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: CPU convergence graph for example saphir.

