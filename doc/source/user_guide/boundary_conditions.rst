.. _boundary_conditions:

-------------------------------
Boundary Conditions
-------------------------------

For all boundaries, it is necessary to provide NekRS with the boundary condition type and, for some types, the boundary condition value.

Boundary condition types are assigned in the ``.par`` file. For each specific field such as velocity or passive scalars (except pressure), the parameter
``boundaryTypeMap`` field setting is typically used (pressure boundary conditions are derived from the velocity boundary conditions).
The boundary conditions are assigned as case-insensitive strings in order of the sideset ID.
By default, nekRS assumes the boundary IDs in your mesh start with 1 and are numbered sequentially.
It is therefore necessary that your mesh has numeric sideset IDs, and slightly preferable that the IDs start from 1 and proceed
in consecutive order. However, NekRS can handle non-consecutive numeric sideset IDs as well.

For example, if your mesh has the following sideset IDs and the corresponding required boundary conditions:

.. csv-table:: Desired Boundary Types
   :align: center
   :header: Boundary ID, Velocity, Temperature

   1,``inlet``,``inlet``
   2,``outlet``,``insulated``
   3,``wall``,``flux``
   4,``sym``,``zerogradient``


you can assign boundary conditions for sidesets in consecutive order (1, 2, 3, 4) as:

.. code-block:: ini

   [VELOCITY]
   boundaryTypeMap = inlet, outlet, wall, sym

   [TEMPERATURE]
   boundaryTypeMap = inlet, insulated, flux, zerogradient

If your boundary IDs are not sequential or do not start with 1, they can be explicitly declared using the ``boundaryIDMap`` key in the ``[MESH]`` block.
The ``boundaryIDMap`` key is a comma-separated list of integers corresponding to the boundary IDs in your mesh.
When using the ``boundaryIDMap`` key, NekRS makes no assumptions regarding the boundary ID values.

.. code-block:: ini

   [MESH]
   boundaryIDMap = 3, 4, 1, 2

   [VELOCITY]
   boundaryTypeMap = wall, sym, inlet, outlet

   [TEMPERATURE]
   boundaryTypeMap = flux, zerogradient, inlet, insulated

Boundary condition values are assigned in the ``.udf`` files (formerly in separate ``.oudf`` files, which are supported for the time-
being but may be deprecated in the future) prior to the user-defined functions (``UDF_Setup`` etc). For the sake of convenience,
the part of the ``.udf`` file where this takes place will be henceforth referred to as the "``oudf`` section". The boundary
condition values are assigned within kernels that execute on the device for every point on the boundary. For common boundary
conditions such as inlet or fixed-pressure outlet boundary conditions, nekRS expects the boundary condition values to be set
in corresponding kernels with pre-defined names such as ``codedFixedValueVelocity`` and ``codedFixedValuePressure`` respectively
(see below for further details). The ``bcData`` structure (``src/nrs/bdry/bcData.h``) provides the members necessary to select
the correct sideset (``bcData->id``), include additional conditions based on geometry (e.g. functions based on x, y, z coordinates
or vector operations performed using normals), and allows users to assign boundary condition values for velocity, scalars,
the mesh solver, or pressure.

Boundary conditions in nekRS are internally interpreted based on certain elementary condition types such as ``zeroValue``,
``codedFixedValue``, ``zeroGradient`` etc. However, depending on the type of field, more convenient and descriptive aliases are also
available. For example, for the velocity field, providing ``inlet`` as a key in the ``boundaryTypeMap`` is acceptable; it
automatically gets parsed as ``codedFixedValue`` when the ``par`` file is read by nekRS. These operations are performed in
``src/core/bdry/bcMap.cpp``.
For pressure, it is possible to set Dirichlet or Neumann conditions at the outlet, if the converse boundary condition
(i.e. Neumann or Dirichlet for velocity, respectively) has been specified for the velocity field. In the absence of an
explicitly-specified boundary condition for pressure, NekRS defaults to a zero gradient boundary condition for
pressure at the outlet. Apart from the outlet, numerically-specified values for pressure boundary conditions are seldom necessary.

The available boundary condition types, along with the identifier strings and aliases, are described in the following sections
for the fluid velocity and pressure, and the temperature and passive scalars.

...........................
Fluid Velocity and Pressure
...........................

Two kinds of boundary conditions are applicable to the fluid velocity: essential (Dirichlet) boundary conditions in which the
velocity is specified, and natural (Neumann) boundary conditions in which the traction is specified. For segments that constitute
the boundary, one of these two types of boundary conditions must be assigned to each component of the fluid velocity.

The fluid boundary condition can be *all Dirichlet* if all velocity components of :math:`{\bf u}` are specified, or it can be
*all Neumann* if all three traction components (:math:`\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}`) are specified on the
boundary face.
where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face.
It can also be *mixed Dirichlet/Neumann* if Dirichlet and Neumann conditions are selected for different velocity components.
Neumann boundary conditions for velocity are assigned differently depending on if the no-stress or full-stress formulation is used.
If the :ref:`no-stress formulation <sec:nostress>` is selected, then traction is not defined on the boundary, rather the individual
components of the velocity gradient are defined.
In this case, any Neumann boundary condition imposed must be homogeneous, i.e. equal to zero, and mixed Dirichlet/Neumann
boundaries must be aligned with one of the Cartesian axes.
These conditions are not required for the full-stress formulation, which assigns the components of
the stress tensor directly.

In general, the boundary condition for pressure satisfies the following equation, unless explicitly specified in the sections below.

 .. math::

  \nabla \cdot \frac{1}{\rho}\nabla p = -\nabla \cdot \frac{D \bf u}{D t} +\nabla \cdot \frac{1}{\rho}\left(\nabla \cdot \boldsymbol{\underline \tau}\right) + \nabla \cdot \bf f

where the stress tensor is given as

 .. math::

   \boldsymbol{\underline \tau} = \mu\left[\nabla {\bf u} + \left(\nabla {\bf u}\right)^T\right]

Inlet (Dirichlet)
`````````````````

Standard Dirichlet boundary condition for velocity. The key values accepted in the ``par`` file for this type of boundary
condition are  ``v``, ``inlet``, or ``codedFixedValue``.

 .. math::

     {\bf u} \cdot {\bf \hat e_x} &= u_x\\
     {\bf u} \cdot {\bf \hat e_y} &= u_y\\
     {\bf u} \cdot {\bf \hat e_z} &= u_z

where, :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face, :math:`{\bf \hat e_x}`, :math:`{\bf \hat e_y}`,
and :math:`{\bf \hat e_z}` are unit vectors aligned with the Cartesian axes and :math:`u_x`, :math:`u_y`, and :math:`u_z` are
set in the ``codedFixedValueVelocity`` kernel within the  ``oudf`` section of the ``.udf`` file, using the ``bc->u/v/w`` members of
the ``bcData`` structure (see the ``turbPipe`` example).

TODO: REMOVE? Inlet (Dirichlet) - local ``vl``
```````````````````````````````````````````````

Standard Dirichlet boundary condition for velocity in local coordinates.

.. math::

     {\bf u} \cdot {\bf \hat e_n} &= u_n\\
     {\bf u} \cdot {\bf \hat e_t} &= u_1\\
     {\bf u} \cdot {\bf \hat e_b} &= u_2

where, :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the
boundary face, and :math:`u_n`, :math:`u_1`, and :math:`u_2` are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.

Wall
````
Dirichlet boundary condition corresponding to a no-slip wall. The key values accepted in the ``par`` file for this type of boundary
condition are ``w``, ``wall``, or ``zerovalue``.

  .. math::

     \bf u = 0

No additional action (say, specifying the ``0`` value in ``oudf``) beyond specifying the boundary value in ``par`` is
necessary for this boundary condition.

Outlet
``````
The open (outflow) boundary condition arises as a natural boundary condition from the variational formulation of Navier-Stokes.
The key values accepted in the ``par`` file for this type of boundary condition are ``o``, ``outlet``, ``outflow``,
and ``zerogradient``.

.. math::

   p = 0

.. csv-table:: 
   :align: center
   :header: no-stress, full-stress
   :widths: 40,40

   :math:`\nabla {\bf u} = 0`,:math:`\boldsymbol{\underline \tau} \cdot {\bf \hat e_n} = 0`

where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face. No specific kernel needs to be called in ``oudf``/``.udf`` for this
boundary condition.

Pressure outlet, specified pressure
```````````````````````````````````

Similar to a standard outlet, but with a specified pressure. This boundary condition is set by specifying a standard outlet
boundary condition for velocity, and using the ``codedFixedValuePressure`` kernel to specify the pressure in ``oudf`` section.

.. math::

   p = p_a

.. csv-table:: 
   :align: center
   :header: no-stress, full-stress
   :widths: 40,40

   :math:`\nabla {\bf u} = 0`,:math:`\boldsymbol{\underline \tau} \cdot {\bf \hat e_n} = 0`

where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face and :math:`p_a` is set in the
``codedFixedValuePressure`` kernel, using the ``bc->p`` member of the ``bcData`` structure (see the ``turbPipe``
example, and the table for other ``bcData`` members below).

Outlet - normal
```````````````

TODO: is ``on`` supported? Comments in ``bcMap.cpp`` indicate otherwise.

Open boundary with zero velocity in the tangent and bitangent directions.
The key values accepted in the ``par`` file for this type of boundary condition are the shorthand forms  ``onx``,
``ony``, or ``onz``; or the expanded forms ``zeroyzvalue/zerogradient``, ``zeroyzvalue/zerogradient``, or
``zeroyzvalue/zerogradient``, corresponding to normals along x, y, or z directions.

.. math::

   p = 0

.. csv-table::
   :align: center
   :header: no-stress, full-stress
   :widths: 40,40

   :math:`\nabla {\bf u}\cdot{\bf \hat e_n} = 0`,:math:`\left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right) \cdot {\bf \hat e_n} = 0`
   :math:`{\bf u} \cdot {\bf \hat e_t} = 0`,:math:`{\bf u} \cdot {\bf \hat e_t} = 0`
   :math:`{\bf u} \cdot {\bf \hat e_b} = 0`,:math:`{\bf u} \cdot {\bf \hat e_b} = 0`

where :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face.
No additional action beyond specifying these parameters in the ``par`` is necessary for this boundary condition.

Pressure outlet - normal
````````````````````````
TODO: is ``on`` supported? Comments in ``bcMap.cpp`` indicate otherwise.

Similar to an outlet - normal boundary, but with a specified pressure.
The key values accepted in the ``par`` file for this type of boundary condition are either ``on`` or ``zerotvalue/zerogradient``.

.. math::

   p = p_a

.. csv-table:: 
   :align: center
   :header: no-stress, full-stress
   :widths: 40,40

   :math:`\nabla {\bf u}\cdot{\bf \hat e_n} = 0`,:math:`\left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right) \cdot {\bf \hat e_n} = 0`
   :math:`{\bf u} \cdot {\bf \hat e_t} = 0`,:math:`{\bf u} \cdot {\bf \hat e_t} = 0`
   :math:`{\bf u} \cdot {\bf \hat e_b} = 0`,:math:`{\bf u} \cdot {\bf \hat e_b} = 0`

where :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors
on the boundary face.
:math:`p_a` is set in the ``oudf`` section within the ``codedFixedValuePressure`` kernel, using the ``bc->p`` member of
the ``bcData`` structure (see below).
If the surface normal vector is not aligned with a principal Cartesian axis, the full-stress formulation must be used by specifying
``equation = navierStokes+variableViscosity`` in the ``[PROBLEMTYPE]`` section of the ``par`` file.

.. _sec:periodicbc:

Periodic
````````

If possible and physically sensible, one can effect great computational efficiency by considering the problem in a single
geometric unit and requiring periodicity of the field variables.

.. math::

   p\left({\bf x}\right) &= p\left({\bf x} + \boldsymbol{\delta}{\bf x}\right)\\
   {\bf u}\left({\bf x}\right) &= {\bf u}\left({\bf x} + \boldsymbol{\delta}{\bf x}\right)

where :math:`\boldsymbol{\delta}{\bf x}` is the offset vector between two periodic faces.
No ``oudf`` kernel is called for this boundary condition type.

Periodic boundaries are a special case where the boundary condition is enforced on the mesh connectivity level. There is no need
to specify periodic boundary conditions in the ``par`` file. Internally, the key set for these boundaries is ``p`` or ``periodic``.
To use periodic boundary conditions, the surface meshes must be conformal.
For third-party meshes they must also have a corresponding pair of boundary ID values which need to be provided
during conversion, i.e. to ``exo2nek`` or ``gmsh2nek``.

TODO: Still accurate? Additionally, the mesh must be at least 3 elements thick in the direction normal to the periodic boundaries.
TODO: check - why is there a key comparison for ``p`` in ``bcMap.cpp`` ? Can you actually set something in ``par``?

See the `meshing section <Meshing>`__ for further details

Symmetry
````````

Symmetric face or a slip wall.
The key values accepted in the ``par`` file for this type of boundary condition are either ``sym``; ``slip``;
or ``zeronvalue/zerogradient`` for boundary face normals in an arbitrary direction, or ``symx/symy/symz``;
``slipx/slipy/slipz``; or  ``zeroxvalue/zerogradient``, ``zeroyvalue/zerogradient``, or ``zerozvalue/zerogradient`` for
faces along the Cartesian axes.

.. math::

   \nabla p \cdot {\bf \hat e_n} = 0

.. csv-table::
   :align: center
   :header: no-stress, full-stress
   :widths: 40,40

   :math:`{\bf u} \cdot {\bf \hat e_n} = 0`,:math:`{\bf u} \cdot {\bf \hat e_n} = 0`
   :math:`\nabla{\bf u}\cdot {\bf \hat e_t} = 0`,:math:`\left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_t} = 0`
   :math:`\nabla{\bf u}\cdot {\bf \hat e_b} = 0`,:math:`\left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_b} = 0`

where :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face.


(TODO: verify) If the surface normal vector is not aligned with a principal Cartesian axis (i.e. a ``sym``, ``slip``, ``zeronvalue/zerogradient`` boundary condition),
the full-stress formulation must be used by specifying ``equation = navierStokes+variableViscosity`` in the ``[PROBLEMTYPE]`` section of the ``par`` file.
No ``oudf`` kernel is called for this boundary condition type.

Traction
````````

This specifies the local neumann boundary conditions for velocity. The shorthand forms specified in the ``par`` file are ``tractionx``, ``tractiony``, ``tractionz``, or ``traction``. These correspond
to ``zeroxvalue/codedfixedgradient``, ``zeroyvalue/codedfixedgradient``, ``zerozvalue/codedfixedgradient`` and ``zeronvalue/codedfixedgradient`` respectively.

Traction boundary conditions of the type ``tractionx/y/z`` specify traction along the Cartesian directions based on the equations:

.. math::

     p &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_x} &= tr_x\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_y} &= tr_y\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_z} &= tr_z

where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face, :math:`{\bf \hat e_x}`, :math:`{\bf \hat e_y}`, and :math:`{\bf \hat e_z}` are unit vectors aligned with the Cartesian axes and :math:`tr_x`, :math:`tr_y`, and :math:`tr_z` are set in the ``oudf`` section within the ``codedFixedGradientVelocity`` kernel using the ``bc->tr1`` (TODO: wild guess, verify) member of the
``bcData`` structure.

The traction can also be specified in local coordinates, using ``traction``. The boundary conditions represented are
  .. math::

     p &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_n} &= tr_n\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_t} &= tr_1\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_b} &= tr_2

where :math:`tr_n`, :math:`tr_1`, and :math:`tr_2` are set in the same ``oudf`` kernel using ``bc->tr1``, ``bc->tr2``, and (TODO: third vector?) members of the ``bcData`` strucutre, and
 :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face available as ``bc->nx/ny/nz``,
``bc->t1x/t1y/t1z``, and ``bc->t2x/t2y/t2z`` respectively.

The full-stress formulation must be used by specifying ``equation = navierStokes+variableViscosity`` in the ``[PROBLEMTYPE]`` section of the ``par`` file for must be used for this boundary type.
See also the ``gabls1`` example in the nekRS repository.

Traction - local, ``sl``
````````````````````````

TODO: does the "Traction" section adequately cover all of this BC?

Similar to traction, but in local coordinates.

  .. math::

     p &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_n} &= tr_n\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_t} &= tr_1\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_b} &= tr_2

where :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face, and :math:`tr_n`, :math:`tr_1`, and :math:`tr_2` are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.
The :ref:`full-stress formulation <sec:fullstress>` must be used for this boundary type.

Traction - horizontal, ``sh``
`````````````````````````````````````
TODO: does the "Traction" section adequately cover all of this BC?

Similar to symmetry, but with specified non-zero traction in the tangent and bitangent directions given in Cartesian coordinates

  .. math::

     {\bf u} \cdot {\bf \hat e_n} &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_x} &= tr_x\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_y} &= tr_y\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_z} &= tr_z

where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face, :math:`{\bf \hat e_x}`, :math:`{\bf \hat e_y}`, and :math:`{\bf \hat e_z}` are unit vectors aligned with the Cartesian axes and :math:`tr_x`, :math:`tr_y`, and :math:`tr_z` are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.
The :ref:`full-stress formulation <sec:fullstress>` must be used for this boundary type.

Traction - horizontal, local, ``shl``
`````````````````````````````````````
TODO: does the "Traction" section adequately cover all of this BC?

Similar to symmetry, but with specified non-zero traction in the tangent and bitangent directions.

  .. math::

     {\bf u} \cdot {\bf \hat e_n} &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_t} &= tr_1\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_b} &= tr_2

where, :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face, and :math:`tr_1` and :math:`tr_2` are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.
The :ref:`full-stress formulation <sec:fullstress>` must be used for this boundary type.

Interpolation Boundaries
````````````````````````

For shared boundaries across Nek-Nek sessions, the ``int`` or ``interpolation`` boundary type is used. See the ``eddyNekNek``
example for further details (TODO: add more detail?)


Moving Mesh Boundary Conditions
```````````````````````````````

For moving mesh problems based on the arbitrary Lagrangian Eulerian framework, such as fluid-structure
interactions or free surface movement, mesh velocity values can be specified for boundaries with the
boundary type of ``mv``, ``codedFixedValue``, or ``codedFixedValue+moving``. The corresponding
``boundaryTypeMap`` in the ``[MESH]`` block should have ``codedfixedvalue``, ``mv``, or ``inlet``
for moving boundaries and ``zerovalue`` for stationary boundaries (TODO: add reference to moving mesh section).

..
TODO: delete Other BCs section once certain nothing here applies
..
Other BCs
..
`````````
..
.. _tab:BCf:

..
.. csv-table:: Other boundary conditions for velocity
   :header: Identifier,Description,Type,Note
   :widths: 5,30,10,55

   ``E`` , "Interior boundary", --, "Denotes faces that connect adjacent elements"
   ``'   '`` , "Empty", --, "Treated as an interior boundary"
   ``int``, "Interpolated (NEKNEK)",       Dirichlet, "Interpolated from the adjacent overset mesh, see: :ref:`neknek`"
   ``p`` , "Periodic", --, "For periodicity within a single element"
   ``mm`` , "Moving mesh",                 --,        "--"
   ``ms`` , "Moving surface",              --,        "--"
   ``msi``, "Moving internal surface",     --,        "--"
   ``mv`` , "Moving boundary",             Dirichlet, "--"
   ``mvn``, "Moving boundary, normal",     Dirichlet, "Zero velocity in non-normal directions"

..
For an axisymmetric flow geometry, the axis boundary condition (``A``) is provided for boundary segments that lie entirely on the axis of symmetry. 
..
This is essentially a symmetry (mixed Dirichlet/Neumann) boundary condition in which the normal velocity and the tangential traction are set to zero.
..
This requires a 2D mesh where the x-axis is the axis of rotation.

..
.. For free-surface boundary segments, the inhomogeneous traction boundary conditions involve both the surface tension coefficient :math:`\sigma` and the mean curvature of the free surface.

..
.. _sec:tempbcs:

...............................
Temperature and Passive Scalars
...............................

The three types of boundary conditions applicable to the temperature are: essential (Dirichlet) boundary condition in which the temperature is specified; natural (Neumann) boundary
condition in which the heat flux is specified; and mixed (Robin) boundary condition in which the heat flux is dependent on the temperature on the boundary.

For segments that constitute a shared boundary between the fluid and the solid domains, one of the above three types of boundary conditions must be assigned to the temperature.

The two types of Robin boundary condition for temperature are: convection boundary conditions for which the heat flux into the domain depends on the heat transfer coefficient :math:`h_{c}` and
the difference between the environmental temperature :math:`T_{\infty}` and the surface temperature; and radiation boundary conditions for which the heat flux into the domain depends on the
Stefan-Boltzmann constant/view-factor product :math:`h_{rad}` and the difference between the fourth power of the environmental temperature :math:`T_{\infty}` and the fourth power of the surface temperature.

The boundary conditions for the passive scalar fields are analogous to those used for the temperature field.
Thus, the temperature boundary conditions and character identifier codes are identical for the passive scalar fields.
The user can specify an independent set of boundary conditions for each passive scalar field.

Specified value (Dirichlet)
```````````````````````````

Standard Dirichlet boundary condition for temperature and passive scalars. Used for inlets, isothermal walls, etc.
Acceptable ``boundaryMapType`` values for this boundary condition are ``t``, ``inlet``, or ``codedfixedvalue``.

.. math::

   T = temp

where :math:`temp` is set in the ``codedFixedValueScalar`` kernel in the ``oudf`` section.

Flux (Neumann), ``f``
`````````````````````

Standard heat flux boundary condition.
Acceptable ``boundaryMapType`` values for this boundary condition are ``f``, ``flux``, or ``codedfixedgradient``.

.. math::

  \lambda\nabla T \cdot {\bf \hat e_n} = flux

where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face and :math:`flux` is set in the ``codedFixedGradientScalar`` kernel in the ``oudf`` section.

Insulated
`````````

Zero-Neumann boundary condition. Used for insulated walls, outlets, symmetry planes, etc.
Acceptable ``boundaryMapType`` values for this boundary condition are ``i``, ``insulated``, ``zeroflux``, or ``zerogradient``.

.. math::

   \lambda \nabla T \cdot {\bf \hat e_n} = 0

where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face.
No ``oudf`` kernel call is required for this boundary condition.

Outflow
```````

Outlet boundary condition but for passive scalars. Mathematically identical to the insulated boundary condition.
Acceptable ``boundaryMapType`` values for this boundary condition are ``o``, ``outflow``, ``outlet``,  or ``zerogradient``.

Interpolation
`````````````

Same as the interpolation boundary condition for velocity but for passive scalars.
Acceptable ``boundaryMapType`` values for this boundary condition are ``int`` or ``interpolation``.

..
  TODO: check if present in nekRS

..
  Newton cooling (convection), ``c``

..
  ``````````````````````````````````

..
  Robin boundary condition for a surface exposed to a fluid at given temperature and heat transfer coefficient.

..
   math::

..
   \lambda \nabla T \cdot {\bf \hat e_n} = h_c\left(T-T_{\infty}\right)

..
  where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face, :math:`h_c` is the convective heat transfer coefficient, and :math:`T_{\infty}` is the ambient temperature.
..
  The convective heat transfer coefficient and ambient temperature are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.

Periodic, ``P``
```````````````

Periodic boundary conditions require that all fields in the simulation are periodic.
This boundary condition works the same way it does for passive scalars, see above.

..
  math::

..
  T \left({\bf x}\right) = T\left({\bf x}+\boldsymbol{\delta}{\bf x}\right)

..
  where :math:`\boldsymbol{\delta}{\bf x}` is the offset vector between two periodic faces.
..
  The ``userbc`` subroutine is not called for this boundary condition type.
..
  See the fluid velocity and pressure :ref:`periodic boundary condition <sec:periodicbc>` for more information.

..
  Radiative cooling, ``r``
..
  ````````````````````````

..
  Robin boundary condition for a surface where radiation heat transfer is significant.

..
  .. math::

..
   \lambda \nabla T \cdot {\bf \hat e_n} = h_{rad}\left(T^4-T_{\infty}^4\right)

..
  where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face, :math:`h_{rad}` is the radiative heat transfer coefficient, and :math:`T_{\infty}` is the ambient temperature.
..
  The radiative heat transfer coefficient and ambient temperature are set in the ``userbc`` subroutine in the :ref:`user file <sec:userbc>`.

..
  Other BCs
..
  `````````

..
  .. _tab:BCt:

..
  .. csv-table:: Other boundary conditions (Temperature and Passive scalars)
   :widths: 5,10,10,75
   :header: Identifier,Description,Type,Note

   ``E``, Interior boundary, --, "--"
   ``'   '`` , "Empty", --, "Treated as an interior boundary"
   ``O``, Outflow, Neumann, "Identical to ``I``"
   ``p``, Periodic, --, "For periodicity within a single element"
   ``SYM``, Symmetry, Neumann, "Identical to ``I``"


............................
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
nekRS always assumes that the density of the two phases are the same (i.e., no Stefan flow); therefore at the melting front, the boundary condition for the fluid velocity is the same as that for a stationary wall, that is, all velocity components are zero.
If no fluid flow is considered, i.e., only the energy equation is solved, then any internal boundary can be defined as a melting front.
The temperature boundary condition at the melting front corresponds to a Dirichlet condition; that is, the entire segment maintains a constant temperature equal to the user-specified melting temperature :math:`T_{melt}` throughout the solution.
In addition, the volumetric latent heat of fusion :math:`\rho L` for the two phases, which is also assumed to be constant, should be specified.

........................
The ``bcData`` structure
........................

The ``bcData`` structure is used in all common boundary condition kernels present in NekRS. It is defined in
``nrs/bdry/bcData.h``. The default argumenti is called ``bc``, and its members are listed below. They are typically
accessed as ``bc->x``, ``bc->id`` etc. They provide useful information for a given Gauss-Lobatto
Legendre point on a boundary. The in most built-in kernels (e.g. ``codedFixedValueVelocity``) is called for
every point on the boundary, so it is important to avoid loops, and to simply use members of this structure
to access point-specific entities.

A list of all available members is given below. Accessing them gives the user the ability to set boundary condition
values for solution fields, check sideset IDs, specify special boundary condition variables like tractions, and to construct functions
dependent on space, geometry, and time.

 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | Name                                                | Type                         | Description                                      |
 +=====================================================+==============================+==================================================+
 | ``idM``                                             | ``int``                      | mesh ID of a given boundary point                |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``fieldOffset``                                     | ``int``                      | offset for a velocity component or a scalar field|
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``id``                                              | ``int``                      | sideset ID                                       |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``time``                                            | ``double``                   | current simulation time                          |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``x``, ``y``, ``z``                                 | ``dfloat``                   | x/y/z value                                      |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``nx``, ``ny``, ``nz``                              | ``dfloat``                   | normals                                          |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``t1x``, ``t1y``, ``t1z``, ``t2x``, ``t2y``, ``t2z`` | ``dfloat``                  | tangential/bitigential directions                |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``tr1``, ``tr2``                                    | ``dfloat``                   | traction in tangential, bitangential directions  |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``u``, ``v``, ``w``                                 | ``dfloat``                   | velocity                                         |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``p``                                               | ``dfloat``                   | pressure                                         |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``uinterp``, ``vinterp``, ``winterp``               | ``dfloat``                   | interpolated velocity values (Nek Nek)           |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``scalarId``                                        | ``int``                      | ID of scalar (based on ``par`` file              |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``s``                                               | ``dfloat``                   | scalar value                                     |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``flux``                                            | ``dfloat``                   | flux value for scalar                            |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``sinterp``                                         | ``dfloat``                   | interpolated scalar value  (Nek Nek)             |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``meshu``, ``meshv``, ``meshw``                     | ``dfloat``                   | x/y/z mesh velocity components                   |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``trans``, ``diff``                                 | ``dfloat``                   | transport/diffusion coeff.                       |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+
 | ``usrwrk``                                          | ``@globalPtr const dfloat*`` | device analogue of ``nrs->usrwrk``               |
 +-----------------------------------------------------+------------------------------+--------------------------------------------------+

Note that ``fieldOffset`` and ``idM`` can be used together to access elements of vector objects like velocity. For example, if you wish
to access the x, y, and z components of velocity at the current boundary point, that can be accomplished using
``U[bc->idM + 0*bc->fieldOffset]``, ``U[bc->idM + 1*bc->fieldOffset]``, and ``U[bc->idM + 2*bc->fieldOffset]`` respectively,
where ``U`` has been passed as an argument of a custom kernel or is transferred from another array such as ``nrs->usrwrk``, which is
useful for transferring data between the host and the device.



..................
Velocity Recycling
..................
One of the challenges with simulating turbulence with high-fidelity is that the flow at the inlet must have turbulent characteristics
to faithfully reproduce the physics of turbulent flow downstream from the inlet. One of the simplest ways to model such turbulence at
the inlet is to use velocity recycling at the inlet  (`see also <https://www.sciencedirect.com/science/article/pii/S0045793009001601>`__).
The flow at least 3-4 characteristic lengths downstream from the inlet is sampled and interpolated back to the inlet, after rescaling
to ensure mass conservation. As turbulence develops in the simulation through the accumulation of numerical errors, the turbulent flow
downstream is reflected back in the inlet velocity field. Similar rescaling and recycling can be applied for passive scalars as well.
For implementation details, see the ``turbPipe`` example.


............
Dong Outflow
............

Backflow at the outlet in turbulent flows can cause CFD simulations to diverge and/or create aphysical results. The current preferred
outlet boundary condition for preventing this is the `Dong outflow <https://www.sciencedirect.com/science/article/pii/S0021999113008504>`__
boundary condition. It applies a pressure distribution directed out of the fluid domain and normal to the outlet, in order to ensure
all fluid near the outlet leaves the domain cleanly without resulting in any backflow. For implementation details, see the
coded ``codedFixedValuePressure`` kernel in the ``oudf`` section of the ``turbPipe`` example.
