.. _meshing:

Meshing
=======

The first step in setting up a case is to have a mesh of the geometry to spatially
discretize the simulation domain with sufficient fidelity to accurately resolve the
physics of the case based on the end-use application with the correct boundary conditions.

NekRS uses its own binary mesh format, ``.re2``, due to its small size and speed of processing
compared to ASCII-based formats. Simple meshes can be created using the ``genbox`` tool. More
complex spatial domains can be imported in `Gmsh's <https://gmsh.info/>`__ ``.msh``, or in the
Exodus II (``.exo``) mesh format. Nek5000's `gmsh2nek <https://github.com/Nek5000/Nek5000/blob/master/tools/gmsh2nek/README.md>`__
and `exo2nek <https://github.com/Nek5000/Nek5000/blob/master/tools/exo2nek/README.md>`__.

This part of the documentation covers some simple examples that illustrate how to use each of the
aforementioned tools, along with a basic overview of the `mesh_t` struct that is commonly used
in NekRS to perform operations such as modifying the mesh and setting boundary conditions.

TODO: verify if ``cgns2nek`` is still fit for purpose.

# by following the instructions in the :ref:`scripts` section. We recommend looking
# at the The README's for the `cgns2nek <https://github.com/Nek5000/Nek5000/blob/master/tools/cgns2nek/README.md>`__,
# or
# tools
# if using them, with further information also available in the 
# `Nek5000 documentation <http://nek5000.github.io/NekDoc/tools.html>`__.

.. _cubit_mesh:

Using ``genbox``
-----------------
Simple meshes can be generated using Nek5000's meshing tool ``genbox``.

Below is a simple example adapted from Nek5000's `MHD example <https://github.com/Nek5000/NekExamples/tree/master/mhd>`__.

.. code-block:: none

    gpf.rea
    -3                    (spatial dimensions: <0 for .re2)
    2                   (number of fields)	: v,T
    #
    #    Circular Polarized Flow of Galloway & Proctor (Dynamo)
    #
    #
    #========================================================
    #
    Box 1
    -4  -4  -8                              (nelx,nely,nelz for Box) 
    0.0 1.0 1.0                             (x0, x1 ratio or xe_i) 
    0.0 1.0 1.0                             (y0, y1 ratio or ye_j) 
    0.0 1.0 1.0                             (z0, z1 ratio or ze_k) 
    P  ,P  ,P  ,P  ,P  ,P                   (cbx0,  cbx1,  cby0,  cby1,  cbz0,  cbz1)  Velocity (3 characters)
    P  ,P  ,P  ,P  ,P  ,P                   (cbx0,  cbx1,  cby0,  cby1,  cbz0,  cbz1)  Temperature


- The first line of this file supplies the name of an existing 3D ``.rea`` file that has the appropriate run parameters (viscosity, timestep size, etc.). These parameters can be modified later, but it is important that ``.rea`` be a 3D file.
- The next line indicates the number of dimensions, and the negative sign indicates that ``genbox`` should generate an ``.re2``
file (whereas a positive sign would generate a legacy ``.rea`` file compatible with Nek5000). For NekRS, you will always
run with a ``-3`` since NekRS performs 3D simulations only, and is not compatible with ``.rea`` files.
- The third line indicates the number of fields for this simulation.
- The next set of lines just shows how one can place comments into a ``genbox`` input file.
- The line that starts with "Box" indicates that a new box is starting, and that the following lines describe a typical box input.  Other possible key characters (the first character of Box, "B") are "C" and "M", more on those later.
- The first line after "Box" specifies the number of elements in the
  :math:`x`, :math:`y`, and :math:`z` directions. The negative sign indicates
  that you want ``genbox`` to automatically generate the element distribution
  along each axis, rather than providing it by hand.  (More on this below.)
- The next 3 lines specify the starting point, end point, and the relative size of each element
  along each Cartesian direciton. The third argument for the ratio can be adjusted to create a graded mesh,
  or you can simply specify the exact coordinates for a given Cartesian direction, having adjusted the sign of the
  corresponding term on the line after ``Box1``. For example, for the y direction:

.. code-block:: none

   -4  4 -8                      Nelx  Nely
   0.0   1.0   1.0               x0  x1   ratio
   0.00 0.25 0.4 0.85 1.0        y0  y1 ... y4

- The last two lines specify `boundary conditions <Boundary conditions>`__ conditions on each side
  of the box for each field


Using ``gmsh2nek``
-----------------
Prior to using ``gmsh2nek``, it is recommended you compile it using a script that is already
included within Nek5000. In addition to the compilers necessary to use `Nek5000 <https://nek5000.github.io/NekDoc/quickstart.html>`__,
``gmsh2nek`` requires ``cmake``. Simply switch to ``path/to/Nek5000/tools``, and execute the script
 ``./maketools gmsh2nek``. Once it is compiled, the executable will be available in ``Nek5000/bin``.
If you added the ``Nek5000/bin`` folder to your ``$PATH`` environment variable as recommended in the `Nek5000 Quickstart
Guide <https://nek5000.github.io/NekDoc/quickstart.html>`__, ``gmsh2nek`` can be used as any other utility added to your
shell environment's ``$PATH`` or as terminal commands can be used.

Before converting a Gmsh ``.msh`` file using ``gmsh2nek``, ensure it is saved in the appropriate format using the following
checklist

[] NekRS requires HEX20 elements. Before exporting your ``.msh`` file, create such elements throughout your
mesh by clicking *Mesh->Set Order 2* in the Gmsh GUI, using the command ``SetOrder`` or passing the option
``-order 2`` to Gmsh in the terminal . Refer to the `Gmsh documentation <https://gmsh.info/doc/texinfo/gmsh.html>`__ for further details.

[] Export your mesh as Version 2 ``.msh`` file. While both ASCII and binary files are supported, the latter is recommended
for large meshes. Do not check any boxes in the export menu when using the GUI. Simply select Version 2 ASCII or Binary from
the drop-down menu and proceed with the export.

[] Setting up periodic boundaries requires a few additional steps. See the `section on sideset and periodic BCs` <Sidesets and applying boundary conditions>__`.

[] Note that NekRS does not support 2D simulations. ``gmsh2nek`` will export 2D meshes only because it has that capability for
use with Nek5000, which can perform 2D simulations, but these are not compatible with NekRS.

[] Ensure your sidesets are set up correctly. See the
 `section on sideset and periodic BCs` <Sidesets and applying boundary conditions>__`

``gmsh2nek`` will guide you through the mesh conversion process, step-by-step, with helpful prompts. You can also
merge a solid domain mesh when prompted, provided the mesh is conformal with the fluid domain mesh. This allows for
the setup of conjugate heat transfer cases.

# The appropriate options from Gmsh are shown in :numref:`fig:gmshopts`.
# 
# .. _fig:gmshopts:
# 
# .. figure:: gmsh2nek/gmshopts.png
#    :align: center
#    :figclass: align-center
# 
#    Recommended options when saving a Gmsh mesh for compatibility with *Nek5000*
#   
# .. Note::
# 
#   Leave both boxes **unchecked** when exporting the mesh from Gmsh.




Using ``exo2nek``
---------------

Similar to ``gmsh2nek``, the Nek5000 tool ``exo2nek`` tool converts ``.exo`` meshes into the native ``.re2`` format.
It is compiled similary using ``./maketools exo2nek``. All features and restrictions from ``gmsh2nek`` apply here as well.
The element types supported by ``exo2nek`` are HEX20, TET4+WEDGE6, TET4+HEX8+WEDGE6, TET10+WEDGE15. This tool can also
create a conjugate heat transfer mesh using two conformal meshes - one for the solid and one for the fluid domain. For
further information on setting up sidesets and boundary conditions, see the
 `section on sideset and periodic BCs` <Sidesets and applying boundary conditions>__`



Creating a mesh with ?????? - CGNS 
----------------------------------

**TODO** Verify cgns meshes are supported; modify introduction accordingly.

.. _cht_mesh:

Conjugate Heat Transfer
-----------------------
TODO: why is this needed when ``gmsh2nek`` can create CHT meshes?

TET2HEX
-----------------------
TODO: for Kirk
As Nek5000 supports only hexahedral elements, exo2nek includes a feature that automatically converts tetrahedral and prism meshes to pure hexahedral meshes. All tetrahedral elements are converted to 4 hexahedral elements and all wedge elements are converted to 3 hexahedral elements. These conversions are supported for both 1st and 2nd order elements.

- TET4 + WEDGE6 –> HEX8
- TET10 + WEDGE15 –> HEX20


Sidesets and applying boundary conditions
------------------------------------------

Setting up sidesets correctly is important for being able to apply boundary conditions correctly
and to avoid compile-time errors resulting from incorrect sideset definitions. For further information
on the types of boundary conditions available, see the `boundary conditions section <Boundary conditions>`__.

- The sidesets are identified by NekRS on the basis of their numerical ID. ``gmsh2nek`` and ``exo2nek``
will detect any text-based IDs, but those are not used by NekRS internally.

- The numerical IDs must start with 1 and must be in a continuous, increasing order of integers with no gaps.

- The number of sidesets in the fluid or solid domain must match the number of ``boundaryID`` entries in the
``.par`` file for the respective solution field's card (e.g.: velocity and temperature).

- Periodic boundary conditions are supported but for translational periodicity only. In other words, the periodic
sideset pairs must lie along the same normal vector and they must be conformal, otherwise ``exo2nek`` and ``gmsh2nek``
will raise errors. Rotational periodicity is currently not supported, however rotational symmetry boundary conditions
can be used for RANS or laminar LES/DNS cases if appropriate for the simulations' goals.

- Note that the `cbc` array used by Nek5000 for setting boundary conditions within the ``usrdat2`` subroutine
of ``usr`` files is **not** used by NekRS for setting boundary conditions. Anything specified using the legacy Nek5000
approach for setting boundary conditions using the `cbc` array will not impact the BCs applied within the simulation.
The use of the ``par`` file array is recommended for setting BCs. The ``cbc`` array is used mainly for certain internal
Nek5000 routines that are bundled with NekRS (e.g. ``torque_calc`` for drag or torque calculations) or for backwards 
compatiblity with legacy Nek5000 features, such as assigning or modifying sideset IDs for meshes generated through ``genbox``,
which have the appropriate BCs (e.g. ``v  ``, ``P  `` etc) in the ``cbc`` array but do not have a ``boundaryID`` assigned.

- Periodic boundary setup for meshes imported through ``gmsh2nek`` or ``exo2nek`` does not require declaring them
as periodic boundaries in the ``.par`` file, but instead requires changing the periodic sideset pair's boundary IDs
to 0. This is because periodic faces are considered internal faces "connected" to another such face on the corresponding
periodic sideset, and as such are treated the same as all other internal faces - with a boundary ID of zero. Non-zero
IDs are reserved for external boundaries or sidesets with non-periodic boundary conditions. After setting these boundary
IDs to 0, it may be necessary to adjust the ID of other sidesets to keep the numbering consistent with the aforementioned
requirements. This can be accomplished through the following code in the ``usrdat2`` subroutine of the ``usr`` file (TODO: udf?)

.. code-block:: fortran

   subroutine usrdat2
   implicit none
   include 'SIZE'
   include 'TOTAL'
   integer e,f,nfaces

   nfaces = 2*ldim

   do e=1,nelt
   do f=1,2*ndim
      if (boundaryID(ifc,iel).eq. 1) then
        boundaryID(ifc,iel) = 1
      else if (boundaryID(ifc,iel).eq. 2) then ! Periodic sideset 1
        boundaryID(ifc,iel) = 0
      else if (boundaryID(ifc,iel) .eq. 3) then ! Periodic sideset 1
        boundaryID(ifc,iel) = 0
      else if (boundaryID(ifc,iel) .eq. 4) then ! Convert Sideset 4 to Sideset 2 to avoid gaps in numbering
        boundaryID(ifc,iel) = 2
      endif
   enddo
   enddo

   return
   end

The ``mesh_t`` struct
---------------------
TODO: verify wrt v24, add some details about bcData and how BCs are set in udf

This section describes commonly-used variables related to the mesh, which are all stored
on data structures of type ``mesh_t``. For the fluid domain, all mesh information is stored
in the ``nrs->mesh`` object, while for scalars such as temperature, mesh information is stored on the
``nrs->cds->mesh`` object. These meshes differ in cases such as conjugate heat transfer, where the
velocity mesh is distinct from the temperature mesh. NekRS performs domain decomposition to ensure
an even split of the mesh across the number of specified MPI tasks, in order to keep the computational load
across all MPI tasks as even as possible. As a result, the mesh gets split into pieces with approximately
equal degrees of freedom across each MPI task. The ``mesh_t`` members therefore correspond to the local
section of the mesh accessed by the given MPI task. For example, ``Nelements`` corresponds to the local
number of elements allotted to the given MPI task.

To keep the following summary table general, the variable names are referred to simply as living on
the ``mesh`` object, without any differentiation between whether that ``mesh`` object is the object on
``nrs`` or ``nrs->cds``.

================== ============================ ================== =================================================
Variable Name      Size                         Device?            Meaning
================== ============================ ================== =================================================
``comm``           1                                               MPI communicator
``device``         1                                               backend device
``dim``            1                                               spatial dimension of mesh
``elementInfo``    ``Nelements``                                   phase of element (0 = fluid, 1 = solid)
``EToB``           ``Nelements * Nfaces``       :math:`\checkmark` boundary ID for each face
``N``              1                                               polynomial order for each dimension
``NboundaryFaces`` 1                                               *total* number of faces on a boundary (rank sum)
``Nelements``      1                                               number of local elements owned by current process
``Nfaces``         1                                               number of faces per element
``Nfp``            1                                               number of quadrature points per face
``Np``             1                                               number of quadrature points per element
``rank``           1                                               parallel process rank
``size``           1                                               size of MPI communicator
``vmapM``          ``Nelements * Nfaces * Nfp`` :math:`\checkmark` quadrature point index for faces on boundaries
``x``              ``Nelements * Np``           :math:`\checkmark` :math:`x`-coordinates of quadrature points
``y``              ``Nelements * Np``           :math:`\checkmark` :math:`y`-coordinates of quadrature points
``z``              ``Nelements * Np``           :math:`\checkmark` :math:`z`-coordinates of quadrature points
================== ============================ ================== =================================================

Miscellaneous Tips
------------------
TODO:
preconditioner, other settings?
comparison to lower order codes: coarseness, aspect ratio targets
h & p refinement, practical approach
