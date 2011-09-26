Newton's Method and Adaptivity (06-nonlinear)
---------------------------------------------

**Git reference:** Tutorial example `06-nonlinear
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P04-adaptivity/06-nonlinear>`_.

In this example we use the nonlinear model problem from example
`P02-nonlinear/02-newton-spline <http://hpfem.org/hermes/doc/src/hermes2d/P02-nonlinear/03-newton-spline.html>`_,

.. math::

    -\nabla \cdot (\lambda(u)\nabla u) - f(x,y) = 0 \ \ \ \mbox{in } \Omega = (-10,10)^2,

equipped with nonhomogeneous Dirichlet boundary conditions 

.. math::

    u(x, y) = (x+10)(y+10)/100 \ \ \ \mbox{on } \partial \Omega.

This time it will be solved using automatic adaptivity. 

Initial solve on coarse mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We solve the nonlinear problem on the coarse mesh once, to obtain a good starting 
point for the Newton's method on the first reference mesh.
First we calculate an initial coefficient vector $\bfY_0$ by projecting 
the initial condition on the coarse mesh::

    // Project the initial condition on the FE space to obtain initial
    // coefficient vector for the Newton's method.
    info("Projecting initial condition to obtain initial vector on the coarse mesh.");
    scalar* coeff_vec_coarse = new scalar[Space::get_num_dofs(&space)] ;
    InitialSolutionHeatTransfer init_sln(&mesh);
    OGProjection::project_global(&space, &init_sln, coeff_vec_coarse, matrix_solver);

Then we solve the problem on the coarse mesh:

.. sourcecode::
    .

    // Newton's loop on the coarse mesh. This is done to obtain a good
    // starting point for the Newton's method on the reference mesh.
    info("Solving on coarse mesh:");
    bool verbose = true;
    bool jacobian_changed = true;
    if (!hermes2d.solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse,
        jacobian_changed, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

.. latexcode::
    .

    // Newton's loop on the coarse mesh. This is done to obtain a good
    // starting point for the Newton's method on the reference mesh.
    info("Solving on coarse mesh:");
    bool verbose = true;
    bool jacobian_changed = true;
    if (!hermes2d.solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse,
        rhs_coarse, jacobian_changed, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) 
        error("Newton's iteration failed.");

The resulting coefficient vector is translated into a coarse mesh solution::

    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solution(coeff_vec_coarse, &space, &sln);

Calculating initial coefficient vector on the reference mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the adaptivity loop, we project the best approximation we have 
to obtain an initial coefficient vector on the new fine mesh.
This step is quite important since the reference space is large, and the 
quality of the initial coefficient vector matters a lot. In the first 
adaptivity step, we use the coarse mesh solution (that's why we have 
computed it), and in all other steps we use the previous fine mesh 
solution:

.. sourcecode::
   .

    // Calculate initial coefficient vector on the reference mesh.
    scalar* coeff_vec = new scalar[Space::get_num_dofs(ref_space)];
    if (as == 1)
    {
      // In the first step, project the coarse mesh solution.
      info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
      OGProjection::project_global(ref_space, &sln, coeff_vec, matrix_solver);
    }
    else
    {
      // In all other steps, project the previous fine mesh solution.
      info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
      OGProjection::project_global(ref_space, &ref_sln, coeff_vec, matrix_solver);
    }

.. latexcode::
   .

    // Calculate initial coefficient vector on the reference mesh.
    scalar* coeff_vec = new scalar[Space::get_num_dofs(ref_space)];
    if (as == 1)
    {
      // In the first step, project the coarse mesh solution.
      info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
      OGProjection::project_global(ref_space, &sln, coeff_vec, matrix_solver);
    }
    else
    {
      // In all other steps, project the previous fine mesh solution.
      info("Projecting previous fine mesh solution to obtain initial vector on new fine 
      mesh.");
      OGProjection::project_global(ref_space, &ref_sln, coeff_vec, matrix_solver);
    }

Sample results
~~~~~~~~~~~~~~

We performed an experiment where we used on the coarse mesh (a) orthogonal projection of the 
fine mesh solution and (b) we solved the nonlinear problem on the coarse mesh. 
We found that this difference does not affect convergence significantly, as 
illustrated in the following convergence comparisons.

(1) Convergence in the number of DOF (with and without Newton solve on the new coarse mesh):

.. figure:: 06-nonlinear/conv_dof_compar.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: DOF convergence graph for tutorial example 01-newton-adapt.

(2) Convergence in CPU time (with and without Newton solve on coarse mesh):

.. figure:: 06-nonlinear/conv_cpu_compar.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: CPU convergence graph for tutorial example 01-newton-adapt.

