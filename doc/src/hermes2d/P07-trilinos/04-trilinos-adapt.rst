Trilinos - Adapt (04-trilinos-adapt)
------------------------------------

**Git reference:** Example `04-trilinos-adapt
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P07-trilinos/04-trilinos-adapt>`_.

The purpose of this example is to show how to use Trilinos while adapting mesh.
We'll use the NOX solver, either using Newton's method or JFNK, with or without 
preconditioning. 

Model problem
~~~~~~~~~~~~~

The underlying problem is benchmark 
`layer-internal <http://hpfem.org/hermes/doc/src/hermes3d/benchmarks/layer-interior.html>`_.

Initial coefficient vector on fine mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dimension of the finite element space changes after each adaptivity 
step. Therefore also the initial coefficient vector for NOX has to change.
Here, for simplicity, we always start from the zero vector::

    // Initial coefficient vector for the Newton's method.  
    scalar* coeff_vec = new scalar[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(scalar));

Indeed this is inefficient since lots of useful information from the previous 
reference solution is thrown away. Correctly one should project the last 
reference solution on the new reference mesh, herewith calculating a much 
better iniial vector. This is waiting for someone to implement.

Initializing NOX on reference mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In each adaptivity step we initialize a new NOX solver and preconditioner:

.. sourcecode::
    .

    // Initialize NOX solver.
    NoxSolver solver(&dp, message_type, "GMRES", "Newton", ls_tolerance, "", flag_absresid, abs_resid, 
                     flag_relresid, rel_resid, max_iters);

    // Select preconditioner.
    RCP<Precond> pc = rcp(new MlPrecond("sa"));
    if (PRECOND)
    {
      if (JFNK) solver.set_precond(pc);
      else solver.set_precond("ML");
    }

.. latexcode::
    .

    // Initialize NOX solver.
    NoxSolver solver(&dp, message_type, "GMRES", "Newton", ls_tolerance, "", flag_absresid,
                     abs_resid, flag_relresid, rel_resid, max_iters);

    // Select preconditioner.
    RCP<Precond> pc = rcp(new MlPrecond("sa"));
    if (PRECOND)
    {
      if (JFNK) solver.set_precond(pc);
      else solver.set_precond("ML");
    }

Assembling and solving on reference mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fine mesh problem is solved and the solution coefficient vector converted
into a Solution. Skipping info outputs and visualization this reads:

.. sourcecode::
    .

    info("Assembling by DiscreteProblem, solving by NOX.");
    solver.set_init_sln(coeff_vec);
    if (solver.solve()) {
      Solution::vector_to_solution(solver.get_solution(), ref_space, &ref_sln);
      info("Number of nonlin iterations: %d (norm of residual: %g)", 
        solver.get_num_iters(), solver.get_residual());
      info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
        solver.get_num_lin_iters(), solver.get_achieved_tol());
    }
    else
      error("NOX failed.");

.. latexcode::
    .

    info("Assembling by DiscreteProblem, solving by NOX.");
    solver.set_init_sln(coeff_vec);
    if (solver.solve()) {
      Solution::vector_to_solution(solver.get_solution(), ref_space, &ref_sln);
      info("Number of nonlin iterations: %d (norm of residual: %g)", 
        solver.get_num_iters(), solver.get_residual());
      info("Total number of iterations in linsolver: %d (achieved tolerance in the last 
           step: %g)", 
        solver.get_num_lin_iters(), solver.get_achieved_tol());
    }
    else
      error("NOX failed.");

Projecting fine mesh solution on coarse mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step is common to all hp-adaptivity algorithms in Hermes::

    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver);


The rest
~~~~~~~~

Now we have a pair of solutions to guide automatic hp-adaptivity, and 
we proceed as in benchmark "layer-internal".


