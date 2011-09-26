Basic-ie-picard
---------------

**Git reference:** Example `basic-ie-picard <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/examples/richards/basic-ie-picard>`_.

Model problem
~~~~~~~~~~~~~

This example is similar to basic-ie-newton except it uses the 
Picard's method in each time step (not Newton's method).

We assume the time-dependent Richard's equation

.. math::
    :label: richards-basic-ie-picard

       C(h) \frac{\partial h}{\partial t} - \nabla \cdot (K(h) \nabla h) - K'(h) \frac{\partial h}{\partial z}= 0

where $C$ and $K$ are non-differentiable and thus in some cases the Newton's method has convergence problems.
For this reason, the Picard's method is a popular approach. If Picards is used, then Andersson acceleration should be employed.

equipped with a Dirichlet, given by the initial condition.

.. math::

     x*(100. - x)/2.5 * y/100 - 1000. + H\underline{\ }OFFSET

The pressure head 'h' is between -1000 and 0. For convenience, we
increase it by an offset H_OFFSET = 1000. In this way we can start
from a zero coefficient vector.

Weak formulation
~~~~~~~~~~~~~~~~

The corresponding weak formulation reads

.. math::

     \int_{\Omega} C(h) \frac{h^{n+1}_{k+1} - h^{n}}{\triangle t} v d\bm{x} + \int_{\Omega} K(h^{n+1}_{k}) \nabla h^{n+1}_{k+1} \cdot \nabla v d\bm{x} - \int_{\Omega} K'(h^{n+1}_{k}) \frac{\partial h^{n+1}_{k+1}}{\partial z} v d\bm{x} = 0.

Here $n$ refers to time steps and $k$ to Picard's iterations within one time step.

Defining weak forms
~~~~~~~~~~~~~~~~~~~

The weak formulation is a combination of custom Jacobian and Residual weak forms::

    CustomWeakFormRichardsIEPicard::CustomWeakFormRichardsIEPicard(double time_step, Solution* h_time_prev, Solution* h_iter_prev) : WeakForm(1)
    {
      // Jacobian.
      CustomJacobian* matrix_form = new CustomJacobian(0, 0, time_step);
      matrix_form->ext.push_back(h_time_prev);
      matrix_form->ext.push_back(h_iter_prev);
      add_matrix_form(matrix_form);

      // Residual.
      CustomResidual* vector_form = new CustomResidual(0, time_step);
      vector_form->ext.push_back(h_time_prev);
      vector_form->ext.push_back(h_iter_prev);
      add_vector_form(vector_form);
    }

Sample results
~~~~~~~~~~~~~~

Solution at t = 0.01 s:

.. image:: basic-ie-newton/basic-ie-newton-0-01s.png
   :align: center
   :scale: 40%
   :alt: sample result

Solution at t = 0.03 s:

.. image:: basic-ie-newton/basic-ie-newton-0-03s.png
   :align: center
   :scale: 40%
   :alt: sample result




