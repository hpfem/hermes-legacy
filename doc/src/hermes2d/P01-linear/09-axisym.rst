Axisymmetric Problems (09-axisym)
---------------------------------

**Git reference:** Tutorial example `09-axisym <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P01-linear/09-axisym>`_. 

Model problem
~~~~~~~~~~~~~

We solve stationary heat transfer in a hollow 
cylindrical object shown in the following schematic picture:

.. figure:: 09-axisym/scheme.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Schematic picture.

The symmetry axis of the object is aligned with the y-axis. The
object stands on a hot plate 

.. math::
  
    u = T_{bottom} \ \ \ \mbox{on}\ \Gamma_{bottom}

where $\Gamma_{bottom}$ denotes its bottom face.
On the rest of the boundary we prescribe a radiation (Newton) 
condition 

.. math::

    -\lambda \frac{\partial u}{\partial n} = \alpha (u - T_{ext}).
    
Here $\lambda$ is the 
thermal conductivity of the material, $\alpha$ the heat transfer
coefficient between the object and the air, and $T_{ext}$ the
exterior air temperature.

Using default weak forms in axisymmetric mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All default weak forms provided by Hermes can be used 
for 2D planar problems, 3D problems that are symmetric 
about the x-axis, and 3D problems that are symmetric 
about the y-axis. The mode is set via the optional 
parameter GeomType in the constructor of the default 
form. Thus the user can think in terms of the planar
formulation of the problem, without having to bother
with the axisymmetric forms of the differential 
operators.

The planar form of the stationary heat transfer equation is

.. math::

    -\mbox{div}(\lambda \, \mbox{grad}\, u) = 0.

Hermes provides DefaultJacobianDiffusion and DefaultResidualDiffusion 
for the diffusion operator 

.. math::

    -\mbox{div}(\lambda \, \mbox{grad}\, u).

These forms have the following constructors::

    DefaultJacobianDiffusion(int i = 0, int j = 0, std::string area = HERMES_ANY, 
                             HermesFunction* coeff = HERMES_ONE,
                             SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    DefaultJacobianDiffusion(int i, int j, Hermes::vector<std::string> areas,
                             HermesFunction* coeff = HERMES_ONE,
                             SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

and

::

    DefaultResidualDiffusion(int i = 0, std::string area = HERMES_ANY,
                             HermesFunction* coeff = HERMES_ONE,
                             GeomType gt = HERMES_PLANAR);

    DefaultResidualDiffusion(int i, Hermes::vector<std::string> areas,
                             HermesFunction* coeff = HERMES_ONE,
                             GeomType gt = HERMES_PLANAR);

Custom weak forms
~~~~~~~~~~~~~~~~~

The weak formulation is custom because of the Newton boundary condition. 
The class header reads

.. sourcecode::
    .

    class CustomWeakFormPoissonNewton : public WeakForm
    {
    public:
      CustomWeakFormPoissonNewton(double lambda, double alpha, double T0, std::string bdy_heat_flux);
    };

.. latexcode::
    .

    class CustomWeakFormPoissonNewton : public WeakForm
    {
    public:
      CustomWeakFormPoissonNewton(double lambda, double alpha, double T0,
                                  std::string bdy_heat_flux);
    };


and the constructor has the form

.. sourcecode::
    .

    CustomWeakFormPoissonNewton::CustomWeakFormPoissonNewton(double lambda, double alpha, double T0, 
							     std::string bdy_heat_flux) : WeakForm(1)
    {
      // Jacobian form - volumetric.
      add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, HERMES_ANY, new HermesFunction(lambda),
								HERMES_SYM, HERMES_AXISYM_Y));

      // Jacobian form - surface.
      add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf(0, 0, bdy_heat_flux, new HermesFunction(alpha),
								  HERMES_AXISYM_Y));

      // Residual forms - volumetric.
      add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, HERMES_ANY, new HermesFunction(lambda),
								HERMES_AXISYM_Y));

      // Residual form - surface.
      add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf(0, bdy_heat_flux, new HermesFunction(alpha),
								HERMES_AXISYM_Y));
      add_vector_form_surf(new WeakFormsH1::DefaultVectorFormSurf(0, bdy_heat_flux, new HermesFunction(-alpha * T0),
                                                                  HERMES_AXISYM_Y));
    };

.. latexcode::
    .

    CustomWeakFormPoissonNewton::CustomWeakFormPoissonNewton(double lambda, double alpha,
                                 double T0, std::string bdy_heat_flux) : WeakForm(1)
    {
      // Jacobian form - volumetric.
      add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, HERMES_ANY,
                                       new HermesFunction(lambda), HERMES_SYM, 
                                       HERMES_AXISYM_Y));

      // Jacobian form - surface.
      add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf(0, 0, bdy_heat_flux,
                                            new HermesFunction(alpha), HERMES_AXISYM_Y));

      // Residual forms - volumetric.
      add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, HERMES_ANY,
                                       new HermesFunction(lambda), HERMES_AXISYM_Y));

      // Residual form - surface.
      add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf(0, bdy_heat_flux,
                                            new HermesFunction(alpha), HERMES_AXISYM_Y));
      add_vector_form_surf(new WeakFormsH1::DefaultVectorFormSurf(0, bdy_heat_flux,
                                            new HermesFunction(-alpha * T0),
                                            HERMES_AXISYM_Y));
    };

Sample results
~~~~~~~~~~~~~~

Results for the values $T_{bottom} = 100$, $T_{ext} = 0$, $\lambda = 386$ and $\alpha = 20$ are shown 
below. We start with the stationary temperature distribution:

.. figure:: 09-axisym/solution.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Temperature.

and the following figure shows the temperature gradient:

.. figure:: 09-axisym/gradient.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Temperature gradient.

