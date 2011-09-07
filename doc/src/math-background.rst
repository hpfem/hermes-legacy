Mathematical Background
-----------------------

The main strength of Hermes is a modern adaptive higher-order 
finite element technology that includes

 * Curvilinear elements.
 * Reduced mesh generation needs.
 * Arbitrary-level hanging nodes.
 * Scalar and vector-valued approximations.
 * Advanced nonlinear solver capabilities.
 * Exponential convergence of adaptive *hp*-FEM.
 * Dozens of time-integration methods. 
 * Adaptivity with dynamical meshes for time-dependent problems.
 * Adaptive multimesh *hp*-FEM for multiphysics coupled problems.

Some the above points are discussed in more detail below:

* **Curvilinear elements**: Approximating curved boundaries or material interfaces via small elements with straight edges belongs to history. It is much more efficient to employ curvilinear elements, such as in the following `acoustics problem <http://hpfem.org/hermes/doc/src/hermes2d/examples/acoustics/horn-axisym.html>`_.

.. figure:: hermes2d/img/acoustic.png
   :align: center
   :scale: 70% 
   :figclass: align-center
   :alt: Illustration of curved elements.

* **Reduced mesh generation needs**: Was the previous result obtained with the mesh below? The answer is *yes*. Do we fail to provide support for traditional mesh generators? The answer is *no*.

.. figure:: hermes2d/img/initmesh.png
   :align: center
   :scale: 40% 
   :figclass: align-center
   :alt: Illustration of coarse initial meshes.

.. raw:: latex

   \newpage

* **Arbitrary-level hanging nodes**: Hermes can handle irregular meshes with arbitrary-level hanging nodes. Why? Because it makes computations much faster.

.. figure:: hermes2d/img/intro/ord_2d_c.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Illustration of arbitrary-level hanging nodes.


* **Exponential convergence of adaptive hp-FEM**: Are you skeptical towards adaptive FEM because it makes things slow? Try adaptive *hp*-FEM. A typical comparison of adaptive low-order FEM and *hp*-FEM is shown below.

.. figure:: hermes2d/img/intro/conv_dof.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Typical convergence curves of FEM with linear and quadratic elements and hp-FEM

Same graphs as above but now in terms of CPU time:

.. figure:: hermes2d/img/intro/conv_cpu.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: CPU convergence graph.

* **Dozens of time-integration methods** are readily available, including the most advanced adaptive implicit higher-order methods. Do not underestimate the time discretization error (below on the left), it can be easily orders of magnitude larger than the error in space (below on the right). 

.. figure:: hermes2d/img/time_error.png
   :align: center
   :scale: 70% 
   :figclass: align-center
   :alt: Time error.

* **Multimesh hp-FEM**: Approximating different physical fields on the same mesh belongs to history. For a given solution component, just one finite element mesh can be optimal. Hermes uses an original adaptive multimesh *hp*-FEM technology to discretize any multiphysics problem *on multiple meshes in a monolithic fashion*. No error due to data transfer between various meshes is present. The following figure illustrates this on a coupled problem of heat and moisture transfer in massive concrete walls of a nuclear reactor vessel. 

.. figure:: hermes2d/img/intro/hm-sln-frame.png
   :align: center
   :scale: 50% 
   :figclass: align-left
   :alt: Illustration of multimesh hp-FEM.

.. figure:: hermes2d/img/intro/hm-mesh-frame.png
   :align: center
   :scale: 50% 
   :figclass: align-right
   :alt: Illustration of multimesh hp-FEM.

.. raw:: html

   <hr style="clear: both; visibility: hidden;">

.. raw:: latex

   \newpage

* **Dynamical meshes for time-dependent problems**: In time-dependent problems, different physical fields or solution components can be approximated on individual meshes that evolve in time independently of each other.

.. figure:: hermes2d/img/intro/flame.jpg
   :align: center
   :scale: 70% 
   :figclass: align-center
   :alt: Adaptive hp-FEM with dynamical meshes for a flame propagation problem. 

* **Wide applicability**: Hermes does not employ any error estimate or another technique that would limit its applicability to some particular class of PDE problems. It allows you to tackle an arbitrary PDE or multiphysics PDE system. Visit the `hp-FEM group home page <http://hpfem.org/>`_ and the `gallery <http://hpfem.org/gallery/>`_ to see examples.

.. figure:: hermes2d/img/intro/ns.jpg
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Image of incompressible viscous flow.

