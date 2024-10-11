.. _models_properties:

Models and Source Terms
=======================

The :ref:`User-Defined Host Functions (.udf) <udf_functions>` file in NekRS provides the necessary 
and sufficient interface to load most of the physics models (and postprocessing capabilities).
For instance, the :ref:`RANS <rans_models>` and :ref:`lowMach compressible <low_mach>` models available in nekRS are
loaded by including corresponding header files in ``.udf`` and by calling the appropriate functions from the standard 
functions in ``.udf``. Appropriate boundary conditions for the momentum and scalar transport equations are specified 
in the :ref:`okl block <okl_block>` (or  in the included ``.oudf`` file). Further, any custom source terms that may need to be added
to the momentum or scalar equations are also interfaced through the ``.udf`` file. 
Before proceeding it is, therefore, highly recommended that users familiarize themselves with all components of :ref:`.udf file <udf_functions>`. 

Turbulence models
-----------------

Large Eddy Simulation (LES)
"""""""""""""""""""""""""""

High pass filter relaxation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

RANS models
"""""""""""

.. _ktau_model:

.. Note::
  RANS model requires two passive scalar fields which must be specified in control parameters ``(.par)``
  file. For details on how to setup the ``.par`` file, refer to the section on :ref:`.par file <parameter_file>` and also refer :ref:`RANS Channel tutorial <tutorial_rans>` for specific example of ``.par`` file setup for :term:`RANS` simulation

The essential routines for the :term:`RANS` models in NekRS are available in the namespace in 
``src/nrs/plugins/RANSktau.hpp``. The default RANS model in nekRS is the :math:`k`-:math:`\tau` model [Tombo2024]_.
Details on the formulation of the :math:`k`-:math:`\tau` can be found :ref:`here <rans_models>`.

To use the :term:`RANS` model in nekRS, first add the necessary include file at the top
of your ``.udf`` file:

.. code-block:: cpp

  #include "RANSktau.hpp"

The header file will make the required :term:`RANS` subroutines accessible in the ``.udf`` file 
which add the necessary source terms for the :math:`k` and :math:`\tau` transport equations and 
modify the diffusion operator in the momentum equation.

Further, in the ``UDF_Setup()`` subroutine, add the following code snippet to initialize the 
:term:`RANS` model,

.. code-block:: cpp
  
  void UDF_Setup()
  {
    nrs->userProperties = &uservp;
    nrs->userScalarSource = &userq;

    const auto ktauFieldStart = 1;

    RANSktau::setup(ktauFieldStart);
  }

``RANSktau::`` is the namespace declared in the header file ``RANSktau.hpp`` which contains all required
:term:`RANS` subroutine call definitions.

``ktauFieldStart`` is the index of the scalar field where the turbulent kinetic energy, ``k``, is stored. In the above example, the :term:`TKE` field corresponds to ``SCALAR01`` as
specified in ``.par`` file (see :ref:`tutorial <tutorial_rans>` for details).

.. warning::
  The ``ktauFieldStart`` index must be consistent with the chosen scalar index specified by user in ``.par`` file for :term:`TKE`.

.. note::
  nekRS assumes that the :math:`\tau` field array always follows the TKE scalar field. Thus, in the above example nekRS assumes :math:`\tau` field index is 2.

``nrs->userProperties`` and ``nrs->userScalarSource`` are the pointer variables to internal subroutines in nekRS
which are used to define the user specified transport properties and source terms for the passive scalar equations, respectively.
As in the above code, these are assigned the pointers to ``uservp`` and ``userq`` routines
which must be defined in the ``.udf`` file as follows,

.. code-block:: cpp

  void uservp(double time)
  {
    RANSktau::updateProperties();
  }

  void userq(double time)
  {
    RANSktau::updateSourceTerms();
  }

The ``updateProperties()`` call computes the diffusion coefficients for the momentum and :math:`k`-:math:`\tau`
equations (see :ref:`RANS theory <rans_models>` for details on RANS model equations),

.. math::
  momentum-equation &\rightarrow \mu + \mu_t \\
  k-equation &\rightarrow \Gamma_k = \mu + \frac{\mu_t}{\sigma_k} \\ 
  \tau-equation &\rightarrow \Gamma_\tau = \mu + \frac{\mu_t}{\sigma_\tau}

.. note::
  ``updateProperties()`` also computes the eddy viscosity, :math:`\mu_t`, required in the above diffusion coefficients. If the user desires to extract :math:`\mu_t` array, say for post-processing purpose, it can be accessed as follows in the ``.udf`` file:
 ``auto o_mue_t = RANSktau::o_mue_t();``

while the ``updateSourceTerms()`` call computes all source terms on the right hand side of the :math:`k` and :math:`\tau` transport equations. 

.. math::
  k-equation &\rightarrow P - \rho \beta^* \frac{k}{\tau} \\
  \tau-equation &\rightarrow -\alpha \rho \tau^2 S^2 + \rho \beta - 8 \Gamma_\tau \left( \nabla \tau^{1/2} \cdot \nabla \tau^{1/2} \right) + C_{D_\tau}

Note that the ``uservp`` and ``userq`` routines are called at each time step by the solver. The above calls will, therefore, update the diffusion properties and source terms at each time step for all GLL points.

The final step in the model setup for the :math:`k`-:math:`\tau` :term:`RANS` model is the specification of the boundary conditions for the :math:`k` and :math:`\tau` transport equations. As explained in the :ref:`RANS theory <rans_models>` section, the wall boundary condition for both :math:`k` and :math:`\tau` equations are zero. These must be explicitly assigned in the :ref:`okl block <okl_block>` section of ``.udf`` file,  

.. code-block:: cpp

  #ifdef __okl__

  void codedFixedValueScalar(bcData *bc)
  {
    if(bc->scalarId == 1 || bc->scalarId == 2) bc->s = 0;
  }

.. note::
  For wall resolved :term:`RANS` simulations, the boundary conditions for both :math:`k` and :math:`\tau` transport equations are of Dirichlet type at the wall and equal to zero.

.. warning::
  It is highly recommended to familiarize with :ref:`okl block <okl_block>` for proper boundary specification. The above example assumes that the computational domain has no inlet boundaries. In case there are inlet boundaries present, they will also have Dirichlet type boundaries for the :math:`k` and :math:`\tau` transport equations and it will be necessary to differentiate the value of :math:`k` and :math:`\tau` at the walls (zero) from those at the inlet (problem dependent). This is done using ``bc->id`` identifier in the :term:`okl block`. 
  
Low-Mach Compressible Model
---------------------------

Custom Source Terms
--------------------

Momentum Equation
"""""""""""""""""

Explicit Source Terms
^^^^^^^^^^^^^^^^^^^^^

Implicit Source Terms
^^^^^^^^^^^^^^^^^^^^^

Scalar Equations
""""""""""""""""

Explicit Source Terms
^^^^^^^^^^^^^^^^^^^^^

Implicit Source Terms
^^^^^^^^^^^^^^^^^^^^^


