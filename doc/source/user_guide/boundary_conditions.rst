.. _boundary_conditions:

-------------------------------
Boundary Conditions
-------------------------------

For all boundaries, it is necessary to provide *NekRS* with the boundary condition *type* and, for some types, the boundary condition *value*.

Boundary condition types are assigned in the ``.par`` file. For each specific field (except pressure), the parameter
``boundaryTypeMap`` field setting. The boundary conditions are assigned as case-insensitive strings in order of the sideset ID,
starting from a sideset ID of 1. It is important that your mesh has numeric sideset IDs and that the IDs start from 1 and proceed
in consecutive order.
Boundary condition values are assigned in the ``.udf`` files (formerly in the separate ``.oudf`` files) prior to the user-defined functions (``UDF_Setup`` etc). The ``bcData`` structure (``src/nrs/bdry/bcData.h``) provides the members necessary to select
the correct sideset (``bcData->id``), include additional conditions based on geometry (e.g. x, y, z coordinates or normals),
and allows users to assign boundary condition values for velocity, scalars, the mesh solver, or pressure.

The basic boundary condition types are ``zeroValue``, ``codedFixedValue``, ``zeroGradient``, ``codedFixedGradient``,
``interpolation``, and ``none``. However, depending on the type of field, more convenient and descriptive aliases are also
available. For example, for the velocity field, providing ``inlet`` as a key in the ``boundaryTypeMap`` is acceptable as it
automatically gets parsed as ``codedFixedValue``. These operations are performed in ``src/core/bdry/bcMap.cpp``. For pressure,
it is possible to set Dirichlet or Neumann conditions at the outlet, if the converse boundary condition (i.e. Neumann or Dirichlet
for velocity, respectively) has been specified for the velocity field. In the absence of an explicitly-specified boundary condition
for pressure, NekRS defaults to a zero gradient boundary condition for pressure at the outlet. Apart from the outlet, pressure boundary
conditions are seldom necessary.

The available boundary condition types, along with the identifier strings and aliases, are described in the following sections for the fluid velocity and pressure, and the temperature and passive scalars.

.. _sec:velbcs:

...........................
Fluid Velocity and Pressure
...........................

Two kinds of boundary conditions are applicable to the fluid velocity: essential (Dirichlet) boundary conditions in which the velocity is specified, and natural (Neumann) boundary conditions in which the traction is specified.
For segments that constitute the boundary :math:`\partial \Omega_f`, see :numref:`fig-walls`, one of these two types of boundary conditions must be assigned to each component of the fluid velocity.

The fluid boundary condition can be *all Dirichlet* if all velocity components of :math:`{\bf u}` are specified, or it can be *all Neumann* if all three traction components (:math:`\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}`) are specified on the boundary face. 
Where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face.
It can also be *mixed Dirichlet/Neumann* if Dirichlet and Neumann conditions are selected for different velocity components.
Neumann boundary conditions for velocity are assigned differently depending on if the no-stress or full-stress formulation is used.
If the :ref:`no-stress formulation <sec:nostress>` is selected, then traction is not defined on the boundary, rather the individual components of the velocity gradient are defined.
In this case, any Neumann boundary condition imposed must be homogeneous, i.e. equal to zero, and mixed Dirichlet/Neumann boundaries must be aligned with one of the Cartesian axes.
These conditions are not required for the :ref:`full-stress formulation <sec:fullstress>`, which assigns the components of the stress tensor directly.

In general, the boundary condition for pressure satisfies the following equation, unless explicitly specified in the sections below.

 .. math::

  \nabla \cdot \frac{1}{\rho}\nabla p = -\nabla \cdot \frac{D \bf u}{D t} +\nabla \cdot \frac{1}{\rho}\left(\nabla \cdot \boldsymbol{\underline \tau}\right) + \nabla \cdot \bf f

where the stress tensor is given as

 .. math::

   \boldsymbol{\underline \tau} = \mu\left[\nabla {\bf u} + \left(\nabla {\bf u}\right)^T\right]

Inlet (Dirichlet), ``v``
````````````````````````

Standard Dirichlet boundary condition for velocity.

 .. math::

     {\bf u} \cdot {\bf \hat e_x} &= u_x\\
     {\bf u} \cdot {\bf \hat e_y} &= u_y\\
     {\bf u} \cdot {\bf \hat e_z} &= u_z
    
Where, :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face, :math:`{\bf \hat e_x}`, :math:`{\bf \hat e_y}`, and :math:`{\bf \hat e_z}` are unit vectors aligned with the Cartesian axes and :math:`u_x`, :math:`u_y`, and :math:`u_z` are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.

Inlet (Dirichlet) - local ``vl``
````````````````````````````````

Standard Dirichlet boundary condition for velocity in local coordinates.

.. math::

     {\bf u} \cdot {\bf \hat e_n} &= u_n\\
     {\bf u} \cdot {\bf \hat e_t} &= u_1\\
     {\bf u} \cdot {\bf \hat e_b} &= u_2
    
Where, :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face, and :math:`u_n`, :math:`u_1`, and :math:`u_2` are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.


Outlet, ``O``
`````````````

The open (outflow) boundary condition arises as a natural boundary condition from the variational formulation of Navier Stokes. 

.. math::

   p = 0

.. csv-table:: 
   :align: center
   :header: no-stress, full-stress
   :widths: 40,40

   :math:`\nabla {\bf u} = 0`,:math:`\boldsymbol{\underline \tau} \cdot {\bf \hat e_n} = 0`

Where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face.
The ``userbc`` subroutine is not called for this boundary condition type.

Pressure outlet, ``o``
``````````````````````

Similar to a standard outlet, but with a specified pressure.

.. math::

   p = p_a

.. csv-table:: 
   :align: center
   :header: no-stress, full-stress
   :widths: 40,40

   :math:`\nabla {\bf u} = 0`,:math:`\boldsymbol{\underline \tau} \cdot {\bf \hat e_n} = 0`

Where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face and :math:`p_a` is set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.

Outlet - normal, ``ON``
```````````````````````

Open boundary with zero velocity in the tangent and bitangent directions.

.. math::

   p = 0

.. csv-table:: 
   :align: center
   :header: no-stress, full-stress
   :widths: 40,40

   :math:`\nabla {\bf u}\cdot{\bf \hat e_n} = 0`,:math:`\left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right) \cdot {\bf \hat e_n} = 0`
   :math:`{\bf u} \cdot {\bf \hat e_t} = 0`,:math:`{\bf u} \cdot {\bf \hat e_t} = 0`
   :math:`{\bf u} \cdot {\bf \hat e_b} = 0`,:math:`{\bf u} \cdot {\bf \hat e_b} = 0`

Where :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face.
If the surface normal vector is not aligned with a principal Cartesian axis, the :ref:`full-stress formulation <sec:fullstress>` must be used.
The ``userbc`` subroutine is not called for this boundary condition type.

Pressure outlet - normal, ``on``
````````````````````````````````

Similar to an outlet - normal boundary, but with a specified pressure.

.. math::

   p = p_a

.. csv-table:: 
   :align: center
   :header: no-stress, full-stress
   :widths: 40,40

   :math:`\nabla {\bf u}\cdot{\bf \hat e_n} = 0`,:math:`\left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right) \cdot {\bf \hat e_n} = 0`
   :math:`{\bf u} \cdot {\bf \hat e_t} = 0`,:math:`{\bf u} \cdot {\bf \hat e_t} = 0`
   :math:`{\bf u} \cdot {\bf \hat e_b} = 0`,:math:`{\bf u} \cdot {\bf \hat e_b} = 0`

Where :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face, and :math:`p_a` is set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.
If the surface normal vector is not aligned with a principal Cartesian axis, the :ref:`full-stress formulation <sec:fullstress>` must be used.

.. _sec:periodicbc:

Periodic, ``P``
```````````````

Where possible, one can effect great computational efficiencies by considering the problem in a single geometric unit and requiring periodicity of the field variables. 

.. math::

   p\left({\bf x}\right) &= p\left({\bf x} + \boldsymbol{\delta}{\bf x}\right)\\
   {\bf u}\left({\bf x}\right) &= {\bf u}\left({\bf x} + \boldsymbol{\delta}{\bf x}\right)

Where :math:`\boldsymbol{\delta}{\bf x}` is the offset vector between two periodic faces.
The ``userbc`` subroutine is not called for this boundary condition type.

Periodic boundaries are a special case where the boundary condition is enforced on the mesh connectivity level. 
To use periodic boundary conditions, the surface meshes must be conformal.
For third-party meshes they must also have a corresponding pair of boundary ID values which need to be provided during conversion, i.e. to ``exo2nek``, ``gmsh2nek``, or ``cgns2nek``. 
Additionally, the mesh must be at least 3 elements thick in the direction normal to the periodic boundaries.

Symmetry, ``SYM``
`````````````````

Symmetric face or a slip wall.

.. math::

   \nabla p \cdot {\bf \hat e_n} = 0

.. csv-table::
   :align: center
   :header: no-stress, full-stress
   :widths: 40,40

   :math:`{\bf u} \cdot {\bf \hat e_n} = 0`,:math:`{\bf u} \cdot {\bf \hat e_n} = 0`
   :math:`\nabla{\bf u}\cdot {\bf \hat e_t} = 0`,:math:`\left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_t} = 0`
   :math:`\nabla{\bf u}\cdot {\bf \hat e_b} = 0`,:math:`\left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_b} = 0`

Where :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face.
If the surface normal vector is not aligned with a principal Cartesian axis, the :ref:`full-stress formulation <sec:fullstress>` must be used.
The ``userbc`` subroutine is not called for this boundary condition type.

Traction, ``s``
```````````````

Full Neumann boundary conditions for velocity.

.. math::

     p &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_x} &= tr_x\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_y} &= tr_y\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_z} &= tr_z

Where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face, :math:`{\bf \hat e_x}`, :math:`{\bf \hat e_y}`, and :math:`{\bf \hat e_z}` are unit vectors aligned with the Cartesian axes and :math:`tr_x`, :math:`tr_y`, and :math:`tr_z` are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.
The :ref:`full-stress formulation <sec:fullstress>` must be used for this boundary type.

Traction - local, ``sl``
````````````````````````

Similar to traction, but in local coordinates.

  .. math::

     p &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_n} &= tr_n\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_t} &= tr_1\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_b} &= tr_2

Where :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face, and :math:`tr_n`, :math:`tr_1`, and :math:`tr_2` are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.
The :ref:`full-stress formulation <sec:fullstress>` must be used for this boundary type.

Traction - horizontal, ``sh``
`````````````````````````````````````

Similar to symmetry, but with specified non-zero traction in the tangent and bitangent directions given in Cartesian coordinates

  .. math::

     {\bf u} \cdot {\bf \hat e_n} &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_x} &= tr_x\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_y} &= tr_y\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_z} &= tr_z

Where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face, :math:`{\bf \hat e_x}`, :math:`{\bf \hat e_y}`, and :math:`{\bf \hat e_z}` are unit vectors aligned with the Cartesian axes and :math:`tr_x`, :math:`tr_y`, and :math:`tr_z` are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.
The :ref:`full-stress formulation <sec:fullstress>` must be used for this boundary type.

Traction - horizontal, local, ``shl``
`````````````````````````````````````

Similar to symmetry, but with specified non-zero traction in the tangent and bitangent directions.

  .. math::

     {\bf u} \cdot {\bf \hat e_n} &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_t} &= tr_1\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_b} &= tr_2

Where, :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face, and :math:`tr_1` and :math:`tr_2` are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.
The :ref:`full-stress formulation <sec:fullstress>` must be used for this boundary type.

Wall, ``W``
```````````

Dirichlet boundary condition corresponding to a no-slip wall.

  .. math::

     \bf u = 0

The ``userbc`` subroutine is not called for this boundary condition type.
  
Other BCs
`````````

.. _tab:BCf:

.. csv-table:: Other boundary conditions for velocity
   :header: Identifier,Description,Type,Note
   :widths: 5,30,10,55

   ``A`` , "Axisymmetric boundary", Mixed, "Can only be used on face 1, treated as ``SYM``, see below"
   ``E`` , "Interior boundary", --, "Denotes faces that connect adjacent elements"
   ``'   '`` , "Empty", --, "Treated as an interior boundary"
   ``int``, "Interpolated (NEKNEK)",       Dirichlet, "Interpolated from the adjacent overset mesh, see: :ref:`neknek`"
   ``p`` , "Periodic", --, "For periodicity within a single element"
   ``mm`` , "Moving mesh",                 --,        "--"
   ``ms`` , "Moving surface",              --,        "--"
   ``msi``, "Moving internal surface",     --,        "--"
   ``mv`` , "Moving boundary",             Dirichlet, "--"
   ``mvn``, "Moving boundary, normal",     Dirichlet, "Zero velocity in non-normal directions"

For an axisymmetric flow geometry, the axis boundary condition (``A``) is provided for boundary segments that lie entirely on the axis of symmetry. 
This is essentially a symmetry (mixed Dirichlet/Neumann) boundary condition in which the normal velocity and the tangential traction are set to zero.
This requires a 2D mesh where the x-axis is the axis of rotation.

.. For free-surface boundary segments, the inhomogeneous traction boundary conditions involve both the surface tension coefficient :math:`\sigma` and the mean curvature of the free surface.

.. _sec:tempbcs:

...............................
Temperature and Passive Scalars
...............................

The three types of boundary conditions applicable to the temperature are: essential (Dirichlet) boundary condition in which the temperature is specified; natural (Neumann) boundary condition in which the heat flux is specified; and mixed (Robin) boundary condition in which the heat flux is dependent on the temperature on the boundary.
For segments that constitute the boundary :math:`\partial \Omega_f' \cup \partial \Omega_s'` (refer to Fig. 2.1), one of the above three types of boundary conditions must be assigned to the temperature.

The two types of Robin boundary condition for temperature are: convection boundary conditions for which the heat flux into the domain depends on the heat transfer coefficient :math:`h_{c}` and the difference between the environmental temperature :math:`T_{\infty}` and the surface temperature; and radiation boundary conditions for which the heat flux into the domain depends on the Stefan-Boltzmann constant/view-factor product :math:`h_{rad}` and the difference between the fourth power of the environmental temperature :math:`T_{\infty}` and the fourth power of the surface temperature.

The boundary conditions for the passive scalar fields are analogous to those used for the temperature field.
Thus, the temperature boundary conditions and character identifier codes are identical for the passive scalar fields.
The user can specify an independent set of boundary conditions for each passive scalar field.

Specified value (Dirichlet), ``t``
``````````````````````````````````

Standard Dirichlet boundary condition for temperature and passive scalars. Used for inlets, isothermal walls, etc.

.. math::

   T = temp

Where :math:`temp` is set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.

Flux (Neumann), ``f``
`````````````````````

Standard heat flux boundary condition.

.. math::

  \lambda\nabla T \cdot {\bf \hat e_n} = flux

Where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face and :math:`flux` is set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.

Insulated, ``I``
````````````````

Zero-Neumann boundary condition. Used for insulated walls, outlets, symmetry planes, etc.

.. math::

   \lambda \nabla T \cdot {\bf \hat e_n} = 0

Where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face.
The ``userbc`` subroutine is not called for this boundary condition type.

Newton cooling (convection), ``c``
``````````````````````````````````

Robin boundary condition for a surface exposed to a fluid at given temperature and heat transfer coefficient.

.. math::

   \lambda \nabla T \cdot {\bf \hat e_n} = h_c\left(T-T_{\infty}\right)

Where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face, :math:`h_c` is the convective heat transfer coefficient, and :math:`T_{\infty}` is the ambient temperature.
The convective heat transfer coefficient and ambient temperature are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.

Periodic, ``P``
```````````````

Periodic boundary conditions require that all fields in the simulation are periodic.

.. math::

   T \left({\bf x}\right) = T\left({\bf x}+\boldsymbol{\delta}{\bf x}\right)

Where :math:`\boldsymbol{\delta}{\bf x}` is the offset vector between two periodic faces.
The ``userbc`` subroutine is not called for this boundary condition type.
See the fluid velocity and pressure :ref:`periodic boundary condition <sec:periodicbc>` for more information.

Radiative cooling, ``r``
````````````````````````

Robin boundary condition for a surface where radiation heat transfer is significant.

.. math::

   \lambda \nabla T \cdot {\bf \hat e_n} = h_{rad}\left(T^4-T_{\infty}^4\right)

Where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face, :math:`h_{rad}` is the radiative heat transfer coefficient, and :math:`T_{\infty}` is the ambient temperature.
The radiative heat transfer coefficient and ambient temperature are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.

Other BCs
`````````

.. _tab:BCt:

.. csv-table:: Other boundary conditions (Temperature and Passive scalars)
   :widths: 5,10,10,75
   :header: Identifier,Description,Type,Note

   ``A``, Axisymmetric boundary, --, "treated as ``I``"
   ``E``, Interior boundary, --, "--"
   ``'   '`` , "Empty", --, "Treated as an interior boundary"
   ``int``, "Interpolated (NEKNEK)", "Dirichlet", "Interpolated from the adjacent overset mesh, see: :ref:`neknek`"
   ``O``, Outflow, Neumann, "Identical to ``I``"
   ``p``, Periodic, --, "For periodicity within a single element"
   ``SYM``, Symmetry, Neumann, "Identical to ``I``"
  
.. ............................
  Internal Boundary Conditions
  ............................

  In the spatial discretization, the entire computational domain is subdivided into macro-elements, the boundary segments shared by any two of these macro-elements in :math:`\Omega_f` and :math:`\Omega_s` are denoted as internal boundaries. 
  For fluid flow analysis with a single-fluid system or heat transfer analysis without change-of-phase, internal boundary conditions are irrelevant as the corresponding field variables on these segments are part of the solution. 
  However, for a multi-fluid system and for heat transfer analysis with change-of-phase, special conditions are required at particular internal boundaries, as described in the following.

  For a fluid system composes of multiple immiscible fluids, the boundary (and hence the identity) of each fluid must be tracked, and a jump in the normal traction exists at the fluid-fluid interface if the surface tension coefficient is nonzero.
  For this purpose, the interface between any two fluids of different identity must be defined as a special type of internal boundary, namely, a fluid layer; and the associated surface tension coefficient also needs to be specified.

  In a heat transfer analysis with change-of-phase, NekRS assumes that both phases exist at the start of the solution, and that all solid-liquid interfaces are specified as special internal boundaries, namely, the melting fronts.
  If the fluid flow problem is considered, i.e., the energy equation is solved in conjunction with the momentum and continuity equations, then only the common boundary between the fluid and the solid (i.e., all or portion of :math:`\partial \overline{\Omega}_f'` in :numref:`fig-walls`) can be defined as the melting front.
  In this case, segments on :math:`\partial \overline{\Omega}_f'` that belong to the dynamic melting/freezing interface need to be specified by the user.
  *NekRS* always assumes that the density of the two phases are the same (i.e., no Stefan flow); therefore at the melting front, the boundary condition for the fluid velocity is the same as that for a stationary wall, that is, all velocity components are zero.
  If no fluid flow is considered, i.e., only the energy equation is solved, then any internal boundary can be defined as a melting front.
  The temperature boundary condition at the melting front corresponds to a Dirichlet condition; that is, the entire segment maintains a constant temperature equal to the user-specified melting temperature :math:`T_{melt}` throughout the solution.
  In addition, the volumetric latent heat of fusion :math:`\rho L` for the two phases, which is also assumed to be constant, should be specified.

.. _sec:settingbcs:

..........................................................
Setting Boundary Conditions Types
..........................................................

Assigning boundary condition types in *NekRS* is handled differently depending on if you are using a third-party meshing tool such as *Gmsh*, *ICEM*, *Cubit*, etc. and importing the mesh with ``exo2nek``, ``gmsh2nek``, or ``cgns2nek``, or if you are using a Nek-native tool such as *preNek* or ``genbox`` (see :ref:`tools_genbox`).
In either case, the boundary condition types are set by assigning the corresponding character identifier code in the character boundary condition array, ``cbc``.
The character boundary condition array itself is described :ref:`here <sec:probvars>` and the supported character codes were described in the sections above for :ref:`momentum <sec:velbcs>` and :ref:`temperature and passive scalars <sec:tempbcs>`.
The differences between Nek-native tools and third-party meshing tools are only in how this array gets set.
For Nek-native tools, this array is read directly from the ``.rea`` or ``.re2`` file, which is set based on input provided to the tool itself.
For third-party meshing tools, the boundary *ID* is set in the tool -- e.g. as a *sideset ID* in *ICEM* -- and this information is propagated to the ``.re2`` (mesh) file.
The ``cbc`` array is later filled at runtime based on the boundary IDs.

The recommended method of setting the boundary condition type from the boundary ID is through the ``.par`` file.
This is done through the ``boundaryTypeMap`` key, which is available for the ``VELOCITY``, ``TEMPERATURE``, and ``SCALARXX`` directives.
By default, *NekRS* assumes the boundary IDs are sequential and start from 1.
If this is not the case, the optional ``boundaryIDMap`` key is available for the ``MESH`` directive.
See :ref:`here <case_files_par>` for more information on the ``.par`` file.
A few simple examples of setting the BC types via the ``.par`` file for a mesh with boundary IDs assigned in a third-party mesher are below.

.. warning::

   Setting the boundary condition types in the ``.par`` file is **NOT** supported in V19 or earlier versions. 

In the simplest example, the mesh has 4 boundaries each with a sequentially numbered boundary ID.

.. csv-table:: Desired Boundary Types
   :align: center
   :header: Boundary ID, Velocity, Temperature

   1,``v``,``t``
   2,``O``,``I``
   3,``W``,``f``
   4,``SYM``,``I``

To set the boundary condition types, the ``boundaryTypeMap`` key is used in the ``.par`` file.
The ``boundaryTypeMap`` key is a comma-separated list of the boundary condition types to be assigned to the domain and is avaialble for the velocity, temperature and passive scalar fields.
The character identifiers can always be used for assignment.
Additionally, some of the common boundary types can be assigned using plain-English equivalents in the ``.par`` file only.
For a list of these see :ref:`here <sec:engidentifiers>`.
By default, *NekRS* assumes the boundary IDs in your mesh start with 1 and are numbered sequentially.
Due to the sequential ordering of the boundary IDs in this example, these boundary types can be set using only the ``boundaryTypeMap`` keys in the ``VELOCITY`` and ``TEMPERATURE`` directives:

.. code-block:: ini

   [VELOCITY]
   boundaryTypeMap = v, O, W, SYM

   [TEMPERATURE]
   boundaryTypeMap = t, I, f, I  

If your boundary IDs are not sequential or do not start with 1, they can be explicitly declared using the ``boundaryIDMap`` key in the ``MESH`` directive.
The ``boundaryIDMap`` key is a comma-separated list of integers corresponding to the boundary IDs in your mesh.
When using the ``boundaryIDMap`` key, *NekRS* makes no assumptions regarding the boundary ID values.

.. code-block:: ini

   [MESH]
   boundaryIDMap = 3, 4, 1, 2

   [VELOCITY]
   boundaryTypeMap = W, SYM, v, O  

   [TEMPERATURE]
   boundaryTypeMap = f, I, t, I


## original placeholder text is below


# Boundary conditions
# ===================
# 
# Boundary conditions for each mesh boundary should normally be set in the 
# :ref:`parameter_file`, using the ``boundaryTypeMap`` parameter. This is used 
# within the ``VELOCITY``, ``TEMPERATURE`` or ``SCALARXX`` sections to set the 
# boundary conditions of the respective components of the case.
# 
# Available Types
# ---------------
# 
# The potential values are summarised in the command line help function for 
# nekRS, ``parHelp.txt`` and below.
# 
# .. literalinclude:: ../../parHelp.txt
#    :language: none
#    :lines: 1-3, 149-179
# 
# .. tip::
# 
#    To setup some cases, you may need to use the ``boundaryIDMap`` parameter of 
#    the ``MESH`` section to apply the ``boundaryTypeMap`` options to the correct
#    boundary IDs of the mesh (By default NekRS assumes that the boundaryIDs start
#    at 1). Below is an example where the four boundary conditions are
#    applied to the boundary IDs 389 (``codeFixedValue``), 231 (``zeroGradient``),
#    4 (``zeroValue``) and 23 (``zeroValue``).
# 
#    .. code-block::
# 
#       [VELOCITY]
#       boundaryTypeMap = codedFixedValue, zeroGradient, zeroValue, zeroValue
# 
#       [MESH]
#       boundaryIDMap = 389, 231, 4, 23
# 
# User Defined Value/Gradient
# """""""""""""""""""""""""""
# 
# If a boundary condition requires a value setting rather than just a type, a 
# suitable function will need to be provided within the ``.udf`` file. The name of
# the function used should be in the form of ``codedFixed`` + ``Value/Gradient`` 
# + ``Velocity/Scalar/Mesh``, E.G. ``codedFixedValueVelocity``, 
# ``codedFixedGradientScalar`` and ``codedFixedValueMesh``.
# 
# All of these functions are passed the ``bcData`` struct which has the following
# parameters available:
# 
# +-----------------------------------------------------+------------------------------+------------------------------+
# | Name                                                | Type                         | Description                  |
# +=====================================================+==============================+==============================+
# | ``idM``                                             | ``int``                      |                              |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``fieldOffset``                                     | ``int``                      |                              |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``id``                                              | ``int``                      |                              |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``time``                                            | ``double``                   |                              |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``x``, ``y``, ``z``                                 | ``dfloat``                   |                              |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``nx``, ``ny``, ``nz``                              | ``dfloat``                   | normals                      |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``1x``, ``t1y``, ``t1z``, ``t2x``, ``t2y``, ``t2z`` | ``dfloat``                   | tangential directions        |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``tr1``, ``tr2``                                    | ``dfloat``                   |                              |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``u``, ``v``, ``w``                                 | ``dfloat``                   | velocity                     |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``p``                                               | ``dfloat``                   |                              |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``uinterp``, ``vinterp``, ``winterp``               | ``dfloat``                   | interpolated velocity values |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``scalarId``                                        | ``int``                      |                              |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``s``                                               | ``dfloat``                   |                              |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``flux``                                            | ``dfloat``                   |                              |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``sinterp``                                         | ``dfloat``                   | interpolated scalar value    |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``meshu``, ``meshv``, ``meshw``                     | ``dfloat``                   |                              |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``trans``, ``diff``                                 | ``dfloat``                   | properties                   |
# +-----------------------------------------------------+------------------------------+------------------------------+
# | ``usrwrk``                                          | ``@globalPtr const dfloat*`` |                              |
# +-----------------------------------------------------+------------------------------+------------------------------+
# 
# Internal / Periodic
# """""""""""""""""""
# 
# None is used when a internal boundary condition is required or a periodic 
# boundary condition has been set as part of the mesh and it does not need to be
# considered as part of the standard processing of boundary conditions.
# 
# Perodicity is linked to the mesh connectivity and is handled by the meshing tool.
# 
# .. NekRS supports periodic boundary conditions. To set up a periodic case, first
# .. you need to run ``exo2nek`` to establish the pairings between the periodic sidesets.
# .. All this information will be prompted on the screen by ``exo2nek``;
# .. You will provide the sideset IDs of the periodic boundaries, a search tolerance
# .. for identifying paired sides, and a translation vector that points from one of the
# .. paired sidesets to the other. For example, if you want to have one periodic surface
# .. that is a :math:`z`-plane at :math:`z=-1.0` that is paired to another :math:`z`-plane
# .. at :math:`z=1.0`, the translation vector would be :math:`(0.0, 0.0, 2.0)`.
# 
# .. _periodic_boundary:
# 
# Turbulent Inflow
# ----------------
# 
# Velocity Recycling
# """"""""""""""""""
# 
