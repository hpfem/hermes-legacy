Nonlinear Example (03-nonlinear)
--------------------------------

**Git reference:** Tutorial example `03-nonlinear 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P03-timedep/03-nonlinear>`_.

Model problem
~~~~~~~~~~~~~

Showing a nonlinear time-dependent problem is not really necessary since 
linear and nonlinear problems are treated in the same way. But let us show one anyway.
Our model problem will be a time-dependent version of the example used in part P02:

.. math::

    \frac{\partial u}{\partial t} - \nabla \cdot (\lambda(u)\nabla u) - f(x,y) = 0.

We prescribe nonhomogeneous Dirichlet boundary conditions 

.. math::

    u(x, y) = (x+10)(y+10)/100 \ \ \ \mbox{on } \partial \Omega,

and the same function is used to define the initial condition. The 
problem will be solved in the square $\Omega = (-10,10)^2$ and time interval $(0, T)$.

Weak formulation
~~~~~~~~~~~~~~~~

In order to use arbitrary Runge-Kutta methods, the time derivative is put on the left-hand 
side of the equation and everything else is put on the right:

.. math::

    \frac{\partial u}{\partial t} = \nabla \cdot (\lambda(u)\nabla u) + f(x,y).

Then the time derivative term is skipped and we write the weak form 
for the stationary residual (right-hand side) only:

.. math::

    F(u, v) = \int_{\Omega} - \lambda(u)\nabla u\cdot \nabla v + fv\, \mbox{d}x\mbox{d}y.

Weak forms
~~~~~~~~~~

After inverting the sign of the custom nonlinearity, 

::

    virtual scalar value(double u) const
    {
      return -1 - pow(u, alpha);
    }

the DefaultWeakFormPoisson can be used::

    // Initialize the weak formulation
    CustomNonlinearity lambda(alpha);
    HermesFunction f(heat_src);
    WeakFormsH1::DefaultWeakFormPoisson wf(HERMES_ANY, &lambda, &f);

Time stepping loop
~~~~~~~~~~~~~~~~~~

The time stepping loop is the same as in the previous example::

    // Time stepping loop:
    double current_time = 0; int ts = 1;
    do
    {
      // Perform one Runge-Kutta time step according to the selected Butcher's table.
      info("Runge-Kutta time step (t = %g s, tau = %g s, stages: %d).",
	   current_time, time_step, bt.get_size());
      bool verbose = true;
      Hermes::vector<Solution*> slns_time_prev;
      slns_time_prev.push_back(&sln_time_prev);
      Hermes::vector<Solution*> slns_time_new;
      slns_time_new.push_back(&sln_time_new);

      if (!runge_kutta.rk_time_step(current_time, time_step, slns_time_prev, slns_time_new, 
                                    true, verbose, NEWTON_TOL, NEWTON_MAX_ITER)) {
	error("Runge-Kutta time step failed, try to decrease time step size.");
      }

      // Update time.
      current_time += time_step;

      // Show the new time level solution.
      char title[100];
      sprintf(title, "Solution, t = %g", current_time);
      sview.set_title(title);
      sview.show(&sln_time_new, HERMES_EPS_VERYHIGH);
      oview.show(&space);

      // Copy solution for the new time step.
      sln_time_prev.copy(&sln_time_new);

      // Increase counter of time steps.
      ts++;
    }
    while (current_time < T_FINAL);

