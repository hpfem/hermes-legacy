Time-Integration with Arbitrary Runge-Kutta Methods (02-runge-kutta)
--------------------------------------------------------------------

**Git reference:** Tutorial example `02-runge-kutta <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P03-timedep/02-runge-kutta>`_. 

This example solves the same model problem as example `01-implicit-euler <http://hpfem.org/hermes/doc/src/hermes2d/P03-timedep/01-implicit-euler.html>`_ but it shows how various Runge-Kutta methods can be used for time stepping. Let us begin with a brief introduction 
to the Runge-Kutta methods and Butcher's tables before we explain implementation details.

Runge-Kutta methods and Butcher's tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Runge-Kutta methods are very nicely described on their `Wikipedia page <http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods>`_
and we recommend the reader to have a brief look there before reading further. In particular, read about the representation 
of Runge-Kutta methods via Butcher's tables, and learn how the table looks like if the method is explicit, diagonally-implicit,
and fully implicit. You should understand basic properties of Runge-Kutta methods such as that explicit methods always need
very small time step or they will blow up, that implicit methods can use much larger time steps, and that among implicit methods, 
diagonally-implicit ones are especially desirable because of relatively low computatonal cost.

Butcher's tables currently available in Hermes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a list of predefined Runge-Kutta methods that can be found 
in the file `hermes_common/tables.cpp <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes_common/tables.cpp>`_.
The names of the tables are self-explanatory. The last number is the order of the 
method (a pair of orders for embedded ones). The second-to-last, if provided, is the number of stages.

Explicit methods::

* Explicit_RK_1, 
* Explicit_RK_2, 
* Explicit_RK_3, 
* Explicit_RK_4. 

Implicit methods::

* Implicit_RK_1, 
* Implicit_Crank_Nicolson_2_2, 
* Implicit_SIRK_2_2, 
* Implicit_ESIRK_2_2, 
* Implicit_SDIRK_2_2, 
* Implicit_Lobatto_IIIA_2_2, 
* Implicit_Lobatto_IIIB_2_2, 
* Implicit_Lobatto_IIIC_2_2, 
* Implicit_Lobatto_IIIA_3_4, 
* Implicit_Lobatto_IIIB_3_4, 
* Implicit_Lobatto_IIIC_3_4, 
* Implicit_Radau_IIA_3_5, 
* Implicit_SDIRK_5_4.

Embedded explicit methods::

* Explicit_HEUN_EULER_2_12_embedded, 
* Explicit_BOGACKI_SHAMPINE_4_23_embedded, 
* Explicit_FEHLBERG_6_45_embedded,
* Explicit_CASH_KARP_6_45_embedded, 
* Explicit_DORMAND_PRINCE_7_45_embedded.

Embedded implicit methods::

* Implicit_SDIRK_CASH_3_23_embedded,
* Implicit_SDIRK_BILLINGTON_3_23_embedded,
* Implicit_ESDIRK_TRBDF2_3_23_embedded, 
* Implicit_ESDIRK_TRX2_3_23_embedded,
* Implicit_SDIRK_CASH_5_24_embedded,
* Implicit_SDIRK_CASH_5_34_embedded,
* Implicit_DIRK_FUDZIAH_7_45_embedded. 

Plus, the user is free to define any Butcher's table of his own.

Model problem
~~~~~~~~~~~~~

For the purpose of using Runge-Kutta methods, the equation has to be 
formulated in such a way that the time derivative stands solo on the left-hand side and 
everything else is on the right

.. math::
    :label: eqvit1

       \frac{\partial T}{\partial t} = \frac{\lambda}{c \varrho} \Delta T.

Weak formulation
~~~~~~~~~~~~~~~~

The temporal derivative is skipped, and weak formulation is only done for the right-hand side

.. math::

     F(T) = - \int_{\Omega} \frac{\lambda}{c \varrho} \nabla T\cdot \nabla v
            + \int_{\Gamma_{air}} \frac{\alpha \lambda}{c \varrho} (T_{ext}(t) - T)v.

This is different from example `01-implicit-euler <http://hpfem.org/hermes/doc/src/hermes2d/P03-timedep/01-implicit-euler.html>`_
where the discretization of the time derivative term was hardwired in the weak formulation. 

The function $F$ above is the stationary residual of the equation (i.e., the weak form of the right-hand side).
Since the Runge-Kutta equations are solved using the Newton's method, the reader may want to look at 
the `Newton's method section <http://hpfem.org/hermes/doc/src/hermes2d/P02-nonlinear/newton-intro.html>`_ before
reading further. The weak form for the Jacobian matrix of the stationary residual is

.. math::

     \frac{\partial F_i(Y)}{\partial y_j} = - \int_{\Omega} \frac{\lambda}{c \varrho} \nabla v_j\cdot \nabla v_i 
                  - \int_{\Gamma_{air}} \frac{\alpha \lambda}{c \varrho} v_j v_i.

Defining weak forms
~~~~~~~~~~~~~~~~~~~

The weak forms are very similar to the previous example, except that the terms 
corresponding to the time derivative are missing, and the rest has an opposite sign:

.. sourcecode::
    .

    class CustomWeakFormHeatRK1 : public WeakForm
    {
    public:
      CustomWeakFormHeatRK1(std::string bdy_air, double alpha, double lambda, double heatcap, double rho,
			    double time_step, double* current_time_ptr, double temp_init, double t_final,
			    Solution* prev_time_sln) : WeakForm(1)
      {
	/* Jacobian */
	// Contribution of the time derivative term.
	add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 0, HERMES_ANY, new HermesFunction(1.0 / time_step)));
	// Contribution of the diffusion term.
	add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, HERMES_ANY, new HermesFunction(lambda / (rho * heatcap))));
	// Contribution of the Newton boundary condition.
	add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf(0, 0, bdy_air, new HermesFunction(alpha / (rho * heatcap))));

	/* Residual */
	// Contribution of the time derivative term.
	add_vector_form(new WeakFormsH1::DefaultResidualVol(0, HERMES_ANY, new HermesFunction(1.0 / time_step)));
	// Contribution of the diffusion term.
	add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, HERMES_ANY, new HermesFunction(lambda / (rho * heatcap))));
	CustomVectorFormVol* vec_form_vol = new CustomVectorFormVol(0, time_step);
	vec_form_vol->ext.push_back(prev_time_sln);
	add_vector_form(vec_form_vol);
	// Contribution of the Newton boundary condition.
	add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf(0, bdy_air, new HermesFunction(alpha / (rho * heatcap))));
	// Contribution of the Newton boundary condition.
	add_vector_form_surf(new CustomVectorFormSurf(0, bdy_air, alpha, rho, heatcap,
			     current_time_ptr, temp_init, t_final));
      };

    private:
      // This form is custom since it contains previous time-level solution.
      class CustomVectorFormVol : public WeakForm::VectorFormVol
      {
      public:
	CustomVectorFormVol(int i, double time_step)
	  : WeakForm::VectorFormVol(i), time_step(time_step) 
	{ 
	}

	virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
	{
	  Func<double>* temp_prev_time = ext->fn[0];
	  return -int_u_v<double, scalar>(n, wt, temp_prev_time, v) / time_step;
	}

	virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
	{
	  Func<Ord>* temp_prev_time = ext->fn[0];
	  return -int_u_v<Ord, Ord>(n, wt, temp_prev_time, v) / time_step;

	}

	double time_step;
      };

      // This form is custom since it contains time-dependent exterior temperature.
      class CustomVectorFormSurf : public WeakForm::VectorFormSurf
      {
      public:
	CustomVectorFormSurf(int i, std::string area, double alpha, double rho, double heatcap,
				    double* current_time_ptr, double temp_init, double t_final)
	  : WeakForm::VectorFormSurf(i, area), alpha(alpha), rho(rho), heatcap(heatcap), current_time_ptr(current_time_ptr),
				     temp_init(temp_init), t_final(t_final) 
	{ 
	}

	virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
	{
	    return -alpha / (rho * heatcap) * temp_ext(*current_time_ptr + time_step) * int_v<double>(n, wt, v);
	}

	virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
	{
	    return -alpha / (rho * heatcap) * temp_ext(*current_time_ptr + time_step) * int_v<Ord>(n, wt, v);
	}

	// Time-dependent exterior temperature.
	template<typename Real>
	Real temp_ext(Real t) const 
	{
	  return temp_init + 10. * sin(2*M_PI*t/t_final);
	}

	double alpha, rho, heatcap, *current_time_ptr, temp_init, t_final;
      };
    };

.. latexcode::
    .

    class CustomWeakFormHeatRK1 : public WeakForm
    {
    public:
      CustomWeakFormHeatRK1(std::string bdy_air, double alpha, double lambda, 
                            double heatcap, double rho, double time_step, double* 
                            current_time_ptr, double temp_init, double t_final, 
                            Solution* prev_time_sln): WeakForm(1)
      {
	/* Jacobian */
	// Contribution of the time derivative term.
	add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 0, HERMES_ANY, 
                                           new HermesFunction(1.0 / time_step)));
	// Contribution of the diffusion term.
	add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, HERMES_ANY,
                                      new HermesFunction(lambda / (rho * heatcap))));
	// Contribution of the Newton boundary condition.
	add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf(0, 0, bdy_air,
                                       new HermesFunction(alpha / (rho * heatcap))));

	/* Residual */
	// Contribution of the time derivative term.
	add_vector_form(new WeakFormsH1::DefaultResidualVol(0, HERMES_ANY, 
                                         new HermesFunction(1.0 / time_step)));
	// Contribution of the diffusion term.
	add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, HERMES_ANY, 
                                      new HermesFunction(lambda / (rho * heatcap))));
	CustomVectorFormVol* vec_form_vol = new CustomVectorFormVol(0, time_step);
	vec_form_vol->ext.push_back(prev_time_sln);
	add_vector_form(vec_form_vol);
	// Contribution of the Newton boundary condition.
	add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf(0, bdy_air,
                                        new HermesFunction(alpha / (rho * heatcap))));
	// Contribution of the Newton boundary condition.
	add_vector_form_surf(new CustomVectorFormSurf(0, bdy_air, alpha, rho, heatcap,
			     current_time_ptr, temp_init, t_final));
      };

    private:
      // This form is custom since it contains previous time-level solution.
      class CustomVectorFormVol : public WeakForm::VectorFormVol
      {
      public:
	CustomVectorFormVol(int i, double time_step)
	  : WeakForm::VectorFormVol(i), time_step(time_step) 
	{ 
	}

	virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> 
                             *v, Geom<double> *e, ExtData<scalar> *ext) const 
	{
	  Func<double>* temp_prev_time = ext->fn[0];
	  return -int_u_v<double, scalar>(n, wt, temp_prev_time, v) / time_step;
	}

	virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                        Geom<Ord> *e, ExtData<Ord> *ext) const 
	{
	  Func<Ord>* temp_prev_time = ext->fn[0];
	  return -int_u_v<Ord, Ord>(n, wt, temp_prev_time, v) / time_step;

	}

	double time_step;
      };

      // This form is custom since it contains time-dependent exterior temperature.
      class CustomVectorFormSurf : public WeakForm::VectorFormSurf
      {
      public:
	CustomVectorFormSurf(int i, std::string area, double alpha, double rho, double 
                             heatcap, double* current_time_ptr, double temp_init, 
                             double t_final)
	  : WeakForm::VectorFormSurf(i, area), alpha(alpha), rho(rho), heatcap(heatcap),
                                     current_time_ptr(current_time_ptr),
				     temp_init(temp_init), t_final(t_final) 
	{ 
	}

	virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                             Geom<double> *e, ExtData<scalar> *ext) const 
	{
	    return -alpha / (rho * heatcap) * temp_ext(*current_time_ptr + time_step) * 
                   int_v<double>(n, wt, v);
	}

	virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const 
	{
	    return -alpha / (rho * heatcap) * temp_ext(*current_time_ptr + time_step) *
                   int_v<Ord>(n, wt, v);
	}

	// Time-dependent exterior temperature.
	template<typename Real>
	Real temp_ext(Real t) const 
	{
	  return temp_init + 10. * sin(2*M_PI*t/t_final);
	}

	double alpha, rho, heatcap, *current_time_ptr, temp_init, t_final;
      };
    };

Selecting a Butcher's table
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unless the user wants to define a Butcher's table on his/her own, he/she can select 
a predefined one - for example a second-order diagonally implicit SDIRK-22
method::

    ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

This is followed in main.cpp by creating an instance of the table::

    ButcherTable bt(butcher_table_type);

Initializing Runge-Kutta time stepping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is done by instantiating the RungeKutta class. Passed are
pointers to the discrete problem, Butcher's table, and 
matrix solver::

    // Initialize Runge-Kutta time stepping.
    RungeKutta runge_kutta(&dp, &bt, matrix_solver);

Time-stepping loop
~~~~~~~~~~~~~~~~~~

The time-stepping loop has the form::

    // Time stepping loop:
    int ts = 1;
    do 
    {
      // Perform one Runge-Kutta time step according to the selected Butcher's table.
      info("Runge-Kutta time step (t = %g s, tau = %g s, stages: %d).", 
	   current_time, time_step, bt.get_size());
      bool jacobian_changed = false;
      bool verbose = true;
      if (!runge_kutta.rk_time_step(current_time, time_step, sln_time_prev, 
				    sln_time_new, jacobian_changed, verbose)) {
	error("Runge-Kutta time step failed, try to decrease time step size.");
      }

      // Show the new time level solution.
      char title[100];
      sprintf(title, "Time %3.2f s", current_time);
      Tview.set_title(title);
      Tview.show(sln_time_new);

      // Copy solution for the new time step.
      sln_time_prev->copy(sln_time_new);

      // Increase current time and time step counter.
      current_time += time_step;
      ts++;
    } 
    while (current_time < T_FINAL);

