.. _theory:

Theory
======

This page provides an overview of the governing equations available in nekRS. nekRS includes models for incompressible flow, a partially compressible low-Mach formulation, the Stokes equations, and the :math:`k`-:math:`\tau` :term:`RANS` equations.

.. _ins_model:

Incompressible Flow
-------------------

The governing equations of incompressible flow in dimensional form are

.. math::

    \rho\left(\frac{\partial\mathbf u}{\partial t} +\mathbf u \cdot \nabla \mathbf u\right) = - \nabla p + \nabla \cdot \boldsymbol{\underline\tau} + \rho {\bf f} \,\, \quad \text{  (Momentum)  }

where :math:`\boldsymbol{\underline\tau}=\mu[\nabla \mathbf u+\nabla \mathbf u^{T}]` and :math:`\mathbf f` is a user defined acceleration.

.. math::

    \nabla \cdot \mathbf u =0 \,\, \quad \text{  (Continuity)  }

If the fluid viscosity is constant in the entire domain, the viscous stress tensor can be contracted using the Laplace operator.
Therefore, one may solve the Navier--Stokes equations in either the full-stress formulation

.. _sec:fullstress:

.. math::

   \nabla \cdot \boldsymbol{\underline\tau}=\nabla \cdot \mu[\nabla \mathbf u+\nabla \mathbf u^{T}]

or the no-stress formulation

.. _sec:nostress:

.. math::

   \nabla \cdot \boldsymbol{\underline\tau}=\mu\Delta \mathbf u

- Variable viscosity and RANS models require the full-stress tensor.
- Constant viscosity leads to a simpler stress tensor, which we refer to as the 'no-stress' formulation.

.. _nondimensional_eqs:

Non-Dimensional Navier-Stokes
"""""""""""""""""""""""""""""

Let us introduce the following non-dimensional variables :math:`\mathbf x^\dagger = \frac{\mathbf x}{L}`, :math:`\mathbf u^\dagger = \frac{u}{U}`, :math:`t^\dagger = \frac{tU}{L}`, and :math:`\mathbf f^\dagger =\frac{\mathbf f L}{U^2}`.
Where :math:`L` and :math:`U` are the (constant) characteristic length and velocity scales, respectively.
For the pressure scale we have two options:

- Convective effects are dominant i.e. high velocity flows :math:`p^\dagger = \frac{p}{\rho_0 U^2}`
- Viscous effects are dominant i.e. creeping flows (Stokes flow) :math:`p^\dagger = \frac{p L}{\mu_0 U}`,

where :math:`\rho_0` and :math:`\mu_0` are constant reference values for density and molecular viscosity, respectively.
For highly convective flows we choose the first scaling of the pressure and obtain the non-dimensional Navier-Stokes in the no-stress formulation:

.. math::

    \frac{\partial \mathbf{u^\dagger}}{\partial t^\dagger} + \mathbf{u^\dagger} \cdot \nabla \mathbf{u^\dagger}\ = -\nabla p^\dagger + \frac{1}{Re}\Delta\mathbf u^\dagger + \mathbf f^\dagger.

For the full-stress formulation, we further introduce the dimensionless viscosity, :math:`\mu^\dagger=\frac{\mu}{\mu_0}`, and obtain:

.. math::

    \frac{\partial \mathbf{u^\dagger}}{\partial t^\dagger} + \mathbf{u^\dagger} \cdot \nabla \mathbf{u^\dagger}\ = -\nabla p^\dagger + \frac{1}{Re}\nabla \cdot \left[ \mu^\dagger \left(\nabla\mathbf u^\dagger + \nabla\mathbf u^{\dagger T}\right)\right] + \mathbf f^\dagger,


where :math:`\mathbf f^\dagger` is the dimensionless user defined forcing function.
The non-dimensional number here is the Reynolds number :math:`Re=\frac{\rho_0 U L}{\mu_0}`.

Stokes Flow
-----------

In the case of flows dominated by viscous effects *NekRS* can solve the reduced Stokes equations

.. math::

    \rho\left(\frac{\partial \mathbf u}{\partial t} \right) = - \nabla p + \nabla \cdot \boldsymbol{\underline\tau} + \rho {\bf f} \,\, , \text{in } \Omega_f \text{  (Momentum)  }

where :math:`\boldsymbol{\underline\tau}=\mu[\nabla \mathbf u+\nabla \mathbf u^{T}]` and

.. math::

    \nabla \cdot \mathbf u =0 \,\, , \text{in } \Omega_f  \text{  (Continuity)  }

As described earlier, we can distinguish between the stress and non-stress formulation according to whether the viscosity is variable or not.
The non-dimensional form of these equations can be obtained using the viscous scaling of the pressure.

.. _intro_energy:

Energy Equation
---------------

In addition to the fluid flow, NekRS offers the capability to solve the energy equation, where temperature is treated as a passive scalar for an incompressible flow.

.. math::

    \rho c_{p} \left( \frac{\partial T}{\partial t} + \mathbf u \cdot \nabla T \right) =
       \nabla \cdot (\lambda \nabla T) + \dot{q} \,\,  \text{  (Energy)  } 

where, :math:`\lambda` is the thermal conductivity and :math:`c_p` is the specific heat at constant pressure.

.. _intro_energy_nondim:

Non-Dimensional Energy / Passive Scalar Equation
""""""""""""""""""""""""""""""""""""""""""""""""

A similar non-dimensionalization as for the flow equations using the non-dimensional variables :math:`\mathbf x^\dagger = \frac{\mathbf x}{L}`,  :math:`\mathbf u^\dagger = \frac{u}{U}`, :math:`t^\dagger = \frac{tU}{L}`, :math:`T^\dagger=\frac{T-T_0}{\delta T}`, and :math:`\lambda^\dagger=\frac{\lambda}{\lambda_0}` leads to

.. math::

    \frac{\partial T^\dagger}{\partial t^\dagger} + \mathbf u^\dagger \cdot \nabla T^\dagger =
      \frac{1}{Pe} \nabla \cdot \lambda^\dagger \nabla T^\dagger + q^\dagger \,\,  \text{  (Energy)  } 

where :math:`q^\dagger=\frac{\dot{q} L}{\rho_0 c_{p_0} U \delta T}` is the dimensionless user defined source term.
The non-dimensional number here is the Peclet number, :math:`Pe=\frac{\rho_0 c_{p_0} U L}{\lambda_0}`.

.. _low_mach:

Low-Mach Compressible Flow
--------------------------

The low-Mach compressible equations are derived from the fully compressible Navier-Stokes equations by filtering the acoustic waves, obtained by splitting the pressure into thermodynamic, :math:`p_t`, and hydrodynamic components :math:`p_1`. The resulting low-Mach compressible governing equations, in dimensional form, are (for complete derivation refer [Tombo1997]_ or [Paulucci1982]_)

.. math::

  \nabla \cdot \mathbf{u} &= \beta_T \frac{D T}{D t} - \kappa \frac{D p_{t}}{D t} = Q\\
  \rho \left(\frac{\partial \mathbf{u}}{\partial t} + \mathbf{u} \cdot \nabla \mathbf{u} \right) &= -\nabla p_1 + \nabla \cdot \mu \left(2 \boldsymbol{\underline{S}} - \frac{2}{3} Q \boldsymbol{\underline{I}} \right) + \rho \mathbf{f} \\
  \rho c_p \frac{D T}{D t} &= \nabla \cdot \lambda \nabla T + \dot{q} + \frac{D p_t}{D t}
  
where, :math:`\boldsymbol{\underline{S}} = \frac{1}{2}(\nabla \mathbf{u} + \nabla \mathbf{u}^T)`, :math:`Q` is the divergence, :math:`\boldsymbol{\underline{I}}` is the identity tensor and :math:`\dot{q}` is the volumetric heat source term.
Thermodynamic pressure is the leading order, spatially invariant, term in pressure expansion while hydrodynamic pressure is the first order term. :math:`\beta_T` is the isobaric expansion coefficient and :math:`\kappa` is the isothermal expansion coefficient,

.. math::

  \beta_T &= \frac{1}{\rho} \left.\frac{D \rho}{D t}\right|_p \\
  \kappa &= \frac{1}{\rho} \left.\frac{D \rho}{D t}\right|_T

.. note::

  :math:`D \bullet/ Dt` is the material derivative. Since :math:`p_t` is spatially invariant, the convective component of its material derivative is zero. Therefore, :math:`D p_t/Dt = dp_t/dt`

.. note::

  For an open domain, the thermodynamic pressure is both spatially and temporally constant, i.e. :math:`dp_t/dt = 0`. This further simplifies the above equation system. However, for a closed system, the thermodynamic pressure, although uniform in space, is subject to changing temporally to enforce mass conservation.

The equation system above is not closed and an equation of state (:term:`EOS`) is required to relate the density to the thermodynamic quantities, :math:`\rho = f(p_t,T)`. Further, dynamic viscosity and thermal conductivity also need to be provided by constitutive relations (e.g., Sutherland's law for gases [Sutherland1893]_).

Introducing the non-dimensional variables as follows,

.. math::

  \mathbf{u}^\dagger = \frac{\mathbf{u}}{U}; \,\, T^\dagger = \frac{T}{T_0}; \,\, \vec{x}^\dagger = \frac{\vec{x}}{L};\,\, p_1^\dagger = \frac{p_1}{\rho U^2};\,\, p_t^\dagger = \frac{p_t}{p_0};\,\, t^\dagger = \frac{t U}{L}; \vec{f}^\dagger = \frac{\vec{f}}{f_0} \\
  \rho^\dagger = \frac{\rho}{\rho_0}; \,\, c_p^\dagger = \frac{c_p}{c_{p0}}; \,\, \lambda^\dagger =\frac{\lambda}{\lambda_0}; \,\, \mu^\dagger = \frac{\mu}{\mu_0}; \,\, \beta_T^\dagger = \frac{\beta_T}{\beta_0}; \,\, \kappa^\dagger = \frac{\kappa}{\kappa_0}; \,\, \dot{q}^\dagger = \frac{\dot{q} L}{\rho_0 c_{p0} T_0 U} 

the low-Mach governing equations are obtained as follows. The continuity equation: 

.. math::

  \nabla \cdot \mathbf{u}^\dagger = \beta_0 T_0 \beta_t^\dagger \frac{D T^\dagger}{D t^\dagger} - \kappa_0 p_0 \kappa^\dagger \frac{d p_t^\dagger}{dt^\dagger} = Q^\dagger

mometum equation,

.. math::

  \rho^\dagger \left(\frac{\partial \mathbf{u}^\dagger}{\partial t^\dagger} + \mathbf{u}^\dagger \cdot \nabla \mathbf{u}^\dagger\right) = - \nabla p_1^\dagger + \nabla \cdot \frac{\mu^\dagger}{Re} \left(2 \boldsymbol{\underline{S}}^\dagger - \frac{2}{3} Q^\dagger \boldsymbol{\underline{I}}\right) + \frac{1}{Fr} \rho^\dagger \mathbf{f}^\dagger

and energy equation,

.. math::

  \rho^\dagger c_p^\dagger \frac{D T^\dagger}{D t^\dagger} = \nabla \cdot \frac{\lambda^\dagger}{Re Pr} \nabla T^\dagger + \dot{q}^\dagger + \frac{p_0}{\rho_0 c_{p0} T_0} \frac{d p_t^\dagger}{d t^\dagger}

where :math:`U` and :math:`L` are the characteristic velocity and length scales. :math:`f_0` is reference magnitude of body force.
:math:`p_0` and :math:`T_0` are the reference pressure and temperature, respectively, and :math:`\rho_0, \mu_0, c_{p0}, \lambda_0, \beta_0, \kappa_0` are the corresponding fluid properties (density, dynamic viscosity, specific heat at contant pressure, conductivity, isobaric expansion coefficient and isothermal expansion coefficient, respectively) at reference conditions. 

:math:`Re=\rho_0 U L/\mu_0` is the Reynolds number, :math:`Pr = \mu_0 c_{p0}/\lambda_0` and :math:`Fr=U^2/f_0 L` are the Reynolds number, Prandtl number and Froude number, defined at reference conditions, respectively.
The equations are closed by corresponding EOS in non-dimensional form, :math:`\rho^\dagger = f(p_t^\dagger,T^\dagger)`.
The above equations represent the lowMach equations in the most general form, applicable to real gases.
Depending on the target application and associated assumptions, several simplifications to the equations are possible.
In the subsequent section, we discuss the simplifications corresponding to the most commonly employed assumption, i.e., ideal gas assumption.

Low-Mach Equations with Ideal Gas Assumption
""""""""""""""""""""""""""""""""""""""""""""

The :term:`EOS` for an ideal gas is,

.. math::

  p_t = \rho R T; \,\, c_p-c_v = R \equiv \frac{R}{c_p} = \frac{\gamma - 1}{\gamma}

where :math:`R` is the ideal gas constant, :math:`c_v` is the specific heat at constant volume and :math:`\gamma = c_p/c_v` is the isentropic expansion factor.
In non-dimensional form, considering the properties at reference conditions for non-dimensionalization (i.e., :math:`p_0 = \rho_0 R T_0` and :math:`\frac{R}{c_{p0}}= \frac{\gamma_0-1}{\gamma_0}`), the :term:`EOS` is simply written,

.. math::

  p_t^\dagger = \rho^\dagger T^\dagger

The expansion coefficients, derived from the EOS, in non-dimensional form are,

.. math::

  \beta_T^\dagger = \frac{1}{T^\dagger} \,\, \kappa^\dagger = \frac{1}{p_t^\dagger}

The resulting governing equations for ideal gas assumption, thus, are,

.. math::

  \nabla \cdot \mathbf{u}^\dagger &= \frac{1}{T^\dagger} \frac{D T^\dagger}{D t^\dagger} - \frac{1}{p_t^\dagger} \frac{d p_t^\dagger}{dt^\dagger} = Q^\dagger \\
  \rho^\dagger \left(\frac{\partial \mathbf{u}^\dagger}{\partial t^\dagger} + \mathbf{u}^\dagger \cdot \nabla \mathbf{u}^\dagger\right) &= - \nabla p_1^\dagger + \nabla \cdot \frac{\mu^\dagger}{Re} \left(2 \boldsymbol{\underline{S}}^\dagger - \frac{2}{3} Q^\dagger \boldsymbol{\underline{I}}\right) + \frac{1}{Fr} \rho^\dagger \mathbf{f}^\dagger \\
  \rho^\dagger c_p^\dagger \frac{D T^\dagger}{D t^\dagger} &= \nabla \cdot \frac{\lambda^\dagger}{Re Pr} \nabla T^\dagger + \dot{q}^\dagger + \frac{\gamma_0-1}{\gamma_0} \frac{d p_t^\dagger}{d t^\dagger}
  
.. note::

  For a calorically perfect ideal gas, :math:`c_p` will be constant and non-dimensional :math:`c_p^\dagger = 1`.

.. note::

  Another often used assumption is to consider dynamic viscosity and thermal conductivity independent of temperature (constant). Thus, :math:`\mu^\dagger` and :math:`\lambda^\dagger` will both be unity, further simplifying the above equations.

.. _rans_models:

RANS Models
-----------

For turbulence modeling :term:`nekRS` offers the two-equation :math:`k`-:math:`\tau` :term:`RANS` model [Tombo2024]_ and its :term:`SST` and :term:`DES` variants [Kumar2024]_.
Linear two-equation RANS models rely on the Bousinessq approximation which relates the Reynolds stress tensor to the mean strain rate, :math:`\boldsymbol{\underline {S}}`, linearly through eddy viscosity.
The time-averaged momentum equation is given as,

.. math::

   \rho \left(\frac{\partial \mathbf u}{\partial t} + \mathbf u \cdot \nabla \mathbf u \right) &=
   - \nabla p + \nabla \cdot \left[ (\mu + \mu_t)
   \left( 2 \boldsymbol{\underline S} -
   \frac{2}{3} Q \boldsymbol{\underline I}\right) \right] \\
   \boldsymbol{\underline S} &= \frac{1}{2} \left( \nabla \mathbf u + \nabla\mathbf{u}^T \right)

where :math:`\mu_t` is the turbulent or eddy viscosity and :math:`\boldsymbol{\underline I}` is an identity tensor.
Currently, nekRS only supports incompressible flow where the divergence constraint, :math:`Q`, is zero,

.. math::

	Q = \nabla \cdot \mathbf u = 0

In two-equation models, the description of the local eddy viscosity is given by two additional transported variables, which provide the velocity and length (or time) scale of turbulence. 
The velocity scale is given by turbulent kinetic energy, :math:`k`, while the choice of the second variable, which provides the length or time scale, depends on the specific two-equation model used. In the :math:`k`-:math:`\tau` model, the second transported variable is :math:`\tau`, which is the inverse of the specific dissipation rate :math:`\omega`, and it provides the local time scale of turbulence.

The :math:`k-\tau` model offers certain favorable characteristics over the :math:`k-\omega` model [Wilcox2008]_, including bounded asymptotic behavior of :math:`\tau` and its source terms and favorable near-wall gradients.
These make it especially suited for high-order codes and complex geometries.
It is, therefore, the preferred two-equation RANS model in NekRS.
The :math:`k-\tau` transport equations are,

.. math::

  \rho\left( \frac{\partial k}{\partial t} + \mathbf u \cdot \nabla k\right) & =
  \nabla \cdot (\Gamma_k \nabla k) + P_k - \rho \beta^* \frac{k}{\tau} \\
  \rho\left( \frac{\partial \tau}{\partial t} + \mathbf u \cdot \nabla\tau\right) & =
  \nabla \cdot (\Gamma_\tau \nabla \tau) - \alpha \frac{\tau}{k}P_k + \rho \beta -
  2\frac{\Gamma_\tau}{\tau} (\nabla \tau \cdot \nabla \tau) + C_{D_\tau}

The diffusion terms are given by

.. math::

  \Gamma_k & = \mu + \frac{\mu_t}{\sigma_k} \\
  \Gamma_\tau & = \mu + \frac{\mu_t}{\sigma_\tau}

where, in the :math:`k-\tau` model the eddy viscosity is given by,

.. math::

  \mu_t = \rho k \tau

The production term is given by

.. math::

  P_k = \mu_t\left( \boldsymbol{\underline S : \underline S} \right)

where ":math:`\boldsymbol :`" denotes the double tensor contraction operator.
The final term in the :math:`\tau` equation is the cross-diffusion term, introduced by [Kok2000]_,

.. math::
  :label: ktau_cd

  C_{D_\tau} =(\rho \sigma_d \tau) \text{min}(\nabla k \cdot \nabla \tau,0)

The above term is especially relevant for external flows.
It eliminates non-physical free-stream dependence of the near-wall :math:`\tau` field (see [Tombo2024]_ for details).

All coefficients in the :math:`k-\tau` model are identical to the standard :math:`k-\omega` model [Wilcox2008]_, given as,

.. math::

  \beta = 0.0708; \,\, \beta^*=0.09; \,\, \alpha=0.52; \,\, \sigma_k= \frac{1}{0.6} \,\, \sigma_\tau=2.0; \,\, \sigma_d=\frac{1}{8}

Further theoretical and implementation details on the :math:`k`-:math:`\tau` model can be found in [Tombo2024]_.
