Hollow Conductor
----------------

**Git reference:** Example `hollow conductor <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/examples/thermoelasticity/hollow-conductor>`_.

Problem description
~~~~~~~~~~~~~~~~~~~

The example deals with a massive hollow conductor is heated by induction and 
cooled by water running inside. We will model this problem using linear thermoelasticity 
equations, where the x-displacement, y-displacement, and the temperature will be approximated 
on individual meshes equipped with mutually independent adaptivity mechanisms. 

The computational domain is shown in the following figure and the details of the geometry can be found 
in the corresponding 
`mesh file <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/examples/thermoelasticity/hollow-conductor/domain.mesh>`_.
It is worth mentioning how the circular arcs are defined via NURBS:

::

    curves =
    {
      { 11, 19, 90 },
      { 10, 15, 90 },
      { 16, 6, 90 },
      { 12, 7, 90 }
    }

The triplet on each line consists of two boundary vertex indices and 
the angle of the circular arc.

.. figure:: img-hollow-conductor/domain.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Domain.

For the equations of linear thermoelasticity and the boundary conditions we refer to the 
paper *P. Solin, J. Cerveny, L. Dubcova, D. Andrs: Monolithic Discretization 
of Linear Thermoelasticity Problems via Adaptive Multimesh hp-FEM*,  
`doi.org/10.1016/j.cam.2009.08.092 <http://dx.doi.org/10.1016/j.cam.2009.08.092>`_.

As usual, the multimesh discretization is initialized by creating the master mesh
via copying the xmesh into ymesh and tmesh.

Sample results
~~~~~~~~~~~~~~

Solution (Von Mises stress):

.. figure:: img-hollow-conductor/mises.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Solution.

Solution (temperature):

.. figure:: img-hollow-conductor/temp.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Solution.

Final meshes for $u_1$, $u_2$ and $T$ (h-FEM with linear elements):

.. figure:: img-hollow-conductor/x-mesh-h1.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Solution.

.. figure:: img-hollow-conductor/y-mesh-h1.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Solution.

.. figure:: img-hollow-conductor/t-mesh-h1.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Solution.

Final meshes for $u_1$, $u_2$ and $T$ (h-FEM with quadratic elements):

.. figure:: img-hollow-conductor/x-mesh-h2.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Solution.

.. figure:: img-hollow-conductor/y-mesh-h2.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Solution.

.. figure:: img-hollow-conductor/t-mesh-h2.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Solution.

Final meshes for $u_1$, $u_2$ and $T$ (h-FEM with quadratic elements):

.. figure:: img-hollow-conductor/x-mesh-hp.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Solution.

.. figure:: img-hollow-conductor/y-mesh-hp.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Solution.

.. figure:: img-hollow-conductor/t-mesh-hp.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Solution.

DOF convergence graphs:

.. figure:: img-hollow-conductor/conv_dof.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. figure:: img-hollow-conductor/conv_cpu.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: CPU convergence graph.

Next let us compare multimesh h-FEM with linear elements with the standard (single-mesh)
h-FEM:

.. figure:: img-hollow-conductor/conv_compar_dof.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: DOF convergence graph.

