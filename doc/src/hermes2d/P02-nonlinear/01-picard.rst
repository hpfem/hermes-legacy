Picard's Method (01-picard)
---------------------------

**Git reference:** Tutorial example `01-picard 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P02-nonlinear/01-picard>`_.

Picard's Method
~~~~~~~~~~~~~~~

The Picard's method is a simple approach to the solution of nonlinear problems
where nonlinear products are linearized by moving part of the nonlinearity 
to the previous iteration level. For example, a nonlinear product of the form 
$g(u)u$ would be linearized as $g(u^n) u^{n+1}$. Let us illustrate this on a 
simple model problem.

Model problem
~~~~~~~~~~~~~

We solve a nonlinear equation

.. math::

    -\nabla \cdot (\lambda(u)\nabla u) = f(x,y), \ \ \ u = u_D \ \mbox{on}\ \partial \Omega.

One possible interpretation of this equation is stationary heat transfer where the thermal
conductivity $\lambda$ depends on the temperature $u$, and $f(x,y)$ are heat sources/losses.
Our domain is a square $\Omega = (-10,10)^2$, $f(x,y) = 1$, and the nonlinearity $\lambda$ has the form 

.. math::

    \lambda(u) = 1 + u^\alpha, \ \ \ \alpha = 4.

Recall that $\lambda$ must be entirely positive or entirely negative for the problem to be solvable
according to the theory. The linearized equation has the form 

.. math::

    -\nabla \cdot (\lambda(u^n)\nabla u^{n+1}) = f(x,y), \ \ \ u = u_D \ \mbox{on}\ \partial \Omega.

The Picard's iteration begins from some initial guess $u^0$, in our case a constant 
function, and runs until a convergence criterion is satisfied. Most widely used 
convergence criteria are the relative error between two consecutive iterations, or 
residual of the equation. In this example we will use the former.

Recall that Hermes uses the Newton's method to solve linear problems. Therefore, the 
linearized equation is written as

.. math::

    -\nabla \cdot (\lambda(u^n)\nabla u^{n+1}) - f(x,y) = 0.

The residual weak form reads

.. math::

    \int_{\Omega} \lambda(u^n) \nabla u^{n+1} \cdot \nabla v \, \mbox{d}\bfx 
    - \int_{\Omega}  f(x,y) v \, \mbox{d}\bfx = 0

where $u^n$ is a given function, $u^{n+1}$ the approximate solution, and $v$
a test function. The weak form of the Jacobian is then

.. math::

    \int_{\Omega} \lambda(u^n) \nabla u \cdot \nabla v \, \mbox{d}\bfx

where $u^n$ is a given function, $u$ a basis function, and $v$ a test function. 

Defining custom nonlinearity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The nonlinearity is defined by subclassing the HermesFunction class::

    class CustomNonlinearity : public HermesFunction
    {
    public:
      CustomNonlinearity(double alpha): HermesFunction()
      {
	this->is_const = false;
	this->alpha = alpha;
      }

      virtual scalar value(double u) const
      {
	return 1 + pow(u, alpha);
      }

      virtual Ord value(Ord u) const
      {
	// If alpha is not an integer, then the function
	// is non-polynomial. 
	// NOTE: Setting Ord to 10 is safe but costly,
	// one could save here by looking at special cases 
	// of alpha. 
	return Ord(10);
      }

      protected:
	double alpha;
    };


Defining initial condition
~~~~~~~~~~~~~~~~~~~~~~~~~~

This example uses a constant initial guess::

    // Initialize previous iteration solution for the Picard's method.
    Solution sln_prev_iter(&mesh, INIT_COND_CONST);


Defining weak forms
~~~~~~~~~~~~~~~~~~~

The weak forms are custom because of the external function 
(previous iteration level solution) that needs to be used::

    // NOTE: The linear problem in each step of the Picard's 
    //       method is solved using the Newton's method.

    class CustomWeakFormPicard : public WeakForm
    {
    public:
      CustomWeakFormPicard(Solution* prev_iter_sln, HermesFunction* lambda, HermesFunction* f) 
	: WeakForm(1)
      {
	// Jacobian (custom because of the external function).
	CustomJacobian* matrix_form = new CustomJacobian(0, 0, lambda);
	matrix_form->ext.push_back(prev_iter_sln);
	add_matrix_form(matrix_form);

	// Residual (custom because of the external function).
	CustomResidual* vector_form = new CustomResidual(0, lambda, f);
	vector_form->ext.push_back(prev_iter_sln);
	add_vector_form(vector_form);
      };

    private:
      class CustomJacobian : public WeakForm::MatrixFormVol
      {
      public:
	CustomJacobian(int i, int j, HermesFunction* lambda) : WeakForm::MatrixFormVol(i, j), lambda(lambda)
	{ 
	}

	virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
			     Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
	{
	  scalar result = 0;
	  for (int i = 0; i < n; i++) 
	  {
	    result += wt[i] * lambda->value(ext->fn[0]->val[i]) 
			    * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
	  }
	  return result;
	}

	virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
			Geom<Ord> *e, ExtData<Ord> *ext) const 
	{
	  Ord result = 0;
	  for (int i = 0; i < n; i++) 
	  {
	    result += wt[i] * lambda->value(ext->fn[0]->val[i]) 
			    * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
	  }
	  return result;
	}

	protected:
	  HermesFunction* lambda;
      };

      class CustomResidual : public WeakForm::VectorFormVol
      {
      public:
	CustomResidual(int i, HermesFunction* lambda, HermesFunction* f) 
	  : WeakForm::VectorFormVol(i), lambda(lambda), f(f) 
	{ 
	}

	virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
			     Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
	{
	  scalar result = 0;
	  for (int i = 0; i < n; i++) 
	  {
	    result += wt[i] * lambda->value(ext->fn[0]->val[i]) 
			    * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
	    result += wt[i] * f->value(e->x[i], e->y[i]) * v->val[i];
	  }
	  return result;
	}

	virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
			Geom<Ord> *e, ExtData<Ord> *ext) const 
	{
	  Ord result = 0;
	  for (int i = 0; i < n; i++) 
	  {
	    result += wt[i] * lambda->value(ext->fn[0]->val[i]) * (u_ext[0]->dx[i] 
			    * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
	    result += wt[i] * f->value(e->x[i], e->y[i]) * v->val[i];
	  }
	  return result;
	}

	private:
	  HermesFunction* lambda;
	  HermesFunction* f;
      };
    };

Note that the previous iteration level solution is accessed through ext->fn[0];

Initializing the weak formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The weak formulation is then initialized in the main.cpp file::

    // Initialize the weak formulation.
    CustomNonlinearity lambda(alpha);
    HermesFunction src(-heat_src);
    CustomWeakFormPicard wf(&sln_prev_iter, &lambda, &src);

Picard's iteration loop
~~~~~~~~~~~~~~~~~~~~~~~

The Picard's iteration is performed simply by::

    bool verbose = true;
    hermes2d.solve_picard(&wf, &space, &sln_prev_iter, matrix_solver, PICARD_TOL, 
  	                  PICARD_MAX_ITER, verbose);

Sample results
~~~~~~~~~~~~~~

Approximate solution $u$ for $\alpha = 4$: 

.. image:: 01-picard/solution.png
   :align: center
   :height: 400
   :alt: result for alpha = 4
