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

Using Nek5000 Meshing Tools
---------------------------

------
genbox
------
Simple meshes can be generated using Nek5000's meshing tool ``genbox``.

Below is a simple mesh based on Nek5000's `MHD example <https://github.com/Nek5000/NekExamples/tree/master/mhd>`__.

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
- The next line indicates the number of dimensions, and the negative sign indicates that ``genbox`` should generate an ``.re2`` file (whereas a positive sign would generate a legacy ``.rea`` file compatible with Nek5000). For NekRS, you will always run with a ``-3`` since NekRS performs 3D simulations only, and is not compatible with ``.rea`` files.
- The third line indicates the number of fields for this simulation.
- The next set of lines just shows how one can place comments into a ``genbox`` input file.
- The line that starts with "Box" indicates that a new box is starting, and that the following lines describe a typical box input.  Other possible key characters (the first character of Box, "B") are "C" and "M", more on those later.
- The first line after "Box" specifies the number of elements in the :math:`x`, :math:`y`, and :math:`z` directions. The negative sign indicates that you want ``genbox`` to automatically generate the element distribution along each axis, rather than providing it by hand.  (More on this below.)
- The next 3 lines specify the starting point, end point, and the relative size of each element  along each Cartesian direciton. The third argument for the ratio can be adjusted to create a graded mesh, or you can simply specify the exact coordinates for a given Cartesian direction, having adjusted the sign of the corresponding term on the line after ``Box1``. For example, for the y direction:

.. code-block:: none

   -4  4 -8                      Nelx  Nely
   0.0   1.0   1.0               x0  x1   ratio
   0.00 0.25 0.4 0.85 1.0        y0  y1 ... y4

- The last two lines specify `boundary conditions <Boundary conditions>`__ conditions on each side of the box for each field.

--------
reatore2
--------

The ``reatore2`` tool is a simple mesh conversion utility that will convert a legacy Nek5000 mesh file (``.rea``) into a NekRS-compatible binary mesh file (``.re2``).
The ``.rea`` format is a combined parameter and mesh format.
Running ``reatore2`` will generate both an ``.re2`` file and a new ``.rea`` file.
The new ``.rea`` file will contain the same parameter settings as the previous ``.rea`` file, but will not include any element information.
Instead, the number of elements listed in the new ``.rea`` file will be negative, indicating that the mesh information is contained in a matching ``.re2`` file.
To use ``reatore2``, simply compile it using ``./maketools reatore2`` as discussed for other meshing utilities like ``genbox``. Launch it from the terminal, and simply enter the name of the ``rea`` file, followed by the desired name for the ``re2`` file.


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
checklist:

- NekRS requires HEX20 elements. Before exporting your ``.msh`` file, create such elements throughout your mesh by clicking *Mesh->Set Order 2* in the Gmsh GUI, using the command ``SetOrder`` or passing the option ``-order 2`` to Gmsh in the terminal . Refer to the `Gmsh documentation <https://gmsh.info/doc/texinfo/gmsh.html>`__ for further details.

- Export your mesh as Version 2 ``.msh`` file. While both ASCII and binary files are supported, the latter is recommended for large meshes. Do not check any boxes in the export menu when using the GUI. Simply select Version 2 ASCII or Binary from the drop-down menu and proceed with the export.

- Setting up periodic boundaries requires a few additional steps. See the section on `sideset and boundary condition setup <Sidesets and applying boundary conditions>`__.

- Note that NekRS does not support 2D simulations. ``gmsh2nek`` will export 2D meshes only because it has that capability for use with Nek5000, which can perform 2D simulations, but these are not compatible with NekRS.

- Ensure your sidesets are set up correctly. See the section on `sidesets and boundary condition setup <Sidesets and applying boundary conditions>`__

``gmsh2nek`` will guide you through the mesh conversion process, step-by-step, with helpful prompts. You can also merge a solid domain mesh when prompted, provided the mesh is conformal with the fluid domain mesh. This allows for the setup of conjugate heat transfer cases.


Using ``exo2nek``
---------------

Similar to ``gmsh2nek``, the Nek5000 tool ``exo2nek`` tool converts ``.exo`` meshes into the native ``.re2`` format.
It is compiled similary using ``./maketools exo2nek``. All features and restrictions from ``gmsh2nek`` apply here as well.
The element types supported by ``exo2nek`` are HEX20, TET4+WEDGE6, TET4+HEX8+WEDGE6, TET10+WEDGE15. This tool can also
create a conjugate heat transfer mesh using two conformal meshes - one for the solid and one for the fluid domain. For
further information on setting up sidesets and boundary conditions, see the section on `applying boundary conditions <Sidesets and applying boundary conditions>`__


Conjugate Heat Transfer
-----------------------

NekRS can simulate conjugate heat transfer when provided with a mesh for the solid domain that is conformal to the fluid domain mesh. While Nek5000 required the use of its internal meshing utility ``prenek`` to merge the solid and fluid meshes,
NekRS does not require the use of ``prenek`` (however, meshes generated by ``prenek`` then converted using ``reatore2`` are still compatible with NekRS). Instead, we recommend the use of ``gmsh2nek`` or ``exo2nek`` to merge the solid domain with
the fluid domain when prompted. To distinguish between the mesh element types, reference the ``nrs->mesh->elementInfo[e]`` or ``nrs->mesh->o_elementInfo[e]`` objects, which store 0 for fluid domain elements, and 1 for solid domain elements.

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
on the types of boundary conditions available, see the page that focuses entirely on the `types of boundary conditions
available in NekRS <Boundary conditions>`__.

- The sidesets are identified by NekRS on the basis of their numerical ID. ``gmsh2nek`` and ``exo2nek`` will detect any text-based IDs, but those are not used by NekRS internally.

- The numerical IDs must start with 1 and must be in a continuous, increasing order of integers with no gaps.

- The number of sidesets in the fluid or solid domain must match the number of ``boundaryID`` entries in the ``.par`` file for the respective solution field's card (e.g.: velocity and temperature).

- Periodic boundary conditions are supported but for translational periodicity only. In other words, the periodic sideset pairs must lie along the same normal vector and they must be conformal, otherwise ``exo2nek`` and ``gmsh2nek`` will raise errors. Rotational periodicity is currently not supported, however rotational symmetry boundary conditions can be used for RANS or laminar LES/DNS cases if appropriate for the simulations' goals.

- Note that the `cbc` array used by Nek5000 for setting boundary conditions within the ``usrdat2`` subroutine of ``usr`` files is **not** used by NekRS for setting boundary conditions. Anything specified using the legacy Nek5000 approach for setting boundary conditions using the `cbc` array will not impact the BCs applied within the simulation. The use of the ``par`` file array is recommended for setting BCs. The ``cbc`` array is used mainly for certain internal Nek5000 routines that are bundled with NekRS (e.g. ``torque_calc`` for drag or torque calculations) or for backwards compatiblity with legacy Nek5000 features, such as assigning or modifying sideset IDs for meshes generated through ``genbox``, which have the appropriate BCs (e.g. ``v  ``, ``P  `` etc) in the ``cbc`` array but do not have a ``boundaryID`` assigned.

- Periodic boundary setup for meshes imported through ``gmsh2nek`` or ``exo2nek`` does not require declaring them as periodic boundaries in the ``.par`` file, but instead requires changing the periodic sideset pair's boundary IDs to 0. This is because periodic faces are considered internal faces "connected" to another such face on the corresponding periodic sideset, and as such are treated the same as all other internal faces - with a boundary ID of zero. Non-zero IDs are reserved for external boundaries or sidesets with non-periodic boundary conditions. After setting these boundary IDs to 0, it may be necessary to adjust the ID of other sidesets to keep the numbering consistent with the aforementioned requirements. This can be accomplished through the following code in the ``usrdat2`` subroutine of the ``usr`` file (TODO: udf?)

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

.. table:: Important ``mesh_t`` members
   :name:  mesh_data

================== ============================ ================== =================================================
Variable Name      Size                         Host or Device?           Meaning
================== ============================ ================== =================================================
``comm``           1                            Host               MPI communicator
``device``         1                            Host               backend device
``dim``            1                            Host               spatial dimension of mesh
``elementInfo``    ``Nelements``                Host               phase of element (0 = fluid, 1 = solid)
``EToB``           ``Nelements * Nfaces``       Both               mapping of elements to type of boundary condition
``N``              1                            Host               polynomial order for each dimension
``NboundaryFaces`` 1                            Host               *total* number of faces on a boundary (rank sum)
``Nelements``      1                            Host               number of local elements owned by current process
``Nfaces``         1                            Host               number of faces per element
``Nfp``            1                            Host               number of quadrature points per face
``Np``             1                            Host               number of quadrature points per element
``rank``           1                            Host               parallel process rank
``size``           1                            Host               size of MPI communicator
``Nfields``        1                            Host               Number of fields passed to the PDE solver
``cht``            1                            Host               conjugate heat transfer status (0 = off, 1 = on)
``vmapM``          ``Nelements * Nfaces * Nfp`` Both               quadrature point index for faces on boundaries
``x``              ``Nelements * Np``           Both               :math:`x`-coordinates of physical quadrature points
``y``              ``Nelements * Np``           Both               :math:`y`-coordinates of physical quadrature points
``z``              ``Nelements * Np``           Both               :math:`z`-coordinates of physical quadrature points
``Nvgeo``          ``<chk>``                    <chk>              Volumetric geometric factors
``Nggeo``          ``<chk>``                    <chk>              Second-order volumetric geometric factors
``vertexNodes``    ``<chk>``                    <chk>              Vertex nodes' indices
``edgeNodes``      ``<chk>``                    <chk>              Edge nodes' indices
``edgeNodes``      ``<chk>``                    <chk>              List of element reference interpolation nodes on element faces
``o_LMM``          ``<chk>``                    Device             Lumped mass matrix
``U``              ``Nelements*Np``             Both               Mesh velocity (often used with ALE solver)
``D``              ``Nelements*Np``             Both               1D Differentiation matrix
``o_vgeo``         ``<chk>``                    Device             Volume geometric factors
``o_sgeo``         ``<chk>``                    Device             Surface geometric factors
================== ============================ ================== =================================================

TODO: example of using these, explain what's happening line by line (grep in examples dir)

Mesh modification in NekRS
--------------------------

TODO: mesh modification 3D example w/o usr file (msre?)

Like Nek5000, NekRS also allows for mesh modification during run time. Static, one-time deformations can be performed during the initialization phase using the ``usr`` file. Typically, the best place to perform such mesh
modifications would be the ``usrdat2`` subroutine. One of the most common deformations performed is scaling the entire mesh by a constant factor. For non-dimensionalization of a case, this factor is often the reciprocal of
the hydraulic diameter. From the ``eddyNekNek`` example:

.. code-block:: fortran

      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'
      include 'CASEDATA'

      n = nx1*ny1*nz1*nelt    !  Rescale mesh to [0,2pi]^2

      call cmult(xm1,P_SCALE,n)
      call cmult(ym1,P_SCALE,n)

      do iel=1,nelt
      do ifc=1,2*ndim
         if (cbc(ifc,iel,1) .eq. 'int') then
           boundaryID(ifc,iel) = 1
         endif
      enddo
      enddo

      return
      end

Many more sophisticated transformations are possible as essentially any point on the mesh can be moved in any manner, provided the user ensures that none of the elements ends up with a negative Jacobian. For a more complex mesh modification demonstration,
see the ``periodicHill`` example.

Dynamic mesh movement is possible using the moving mesh solver based on the Arbitrary Lagrangian-Eulerian (ALE) framework. When supplied with a mesh velocity on a moving boundary, the solver automatically moves the mesh according to ALE equations and the
mesh deformation diffusion parameter which controls the spatial blending of the deformation into the entire domain from the moving mesh boundary. For more details, consult the ``mv_cyl`` example. Apart from the ALE solver, the user has the option to control
the mesh deformation themselves using ``solver = user`` in the ``[MESH]`` block, but it becomes the user's responsibility to account for the effects of the mesh movement on the physics of the fluid and any scalars being simulated.

Miscellaneous Tips
------------------
TODO:
preconditioner, other settings?

In higher-order finite elements, mesh refinement can be performed in two ways:

- h-refinement (the addition of elements)
- p-refinement (increasing the polynomial order of all elements).

In practice, we rely on both to achieve optimal mesh resolution. Areas where the flow cross section changes abruptly, where flow separation occurs, or other regions of flow transition such as reactor plena are good candidates for h-refinement. After getting
convergent results with a low polynomial order (say, ``N=3``), it is a good idea to switch to p-refinement to see if the quantity of interest that you have chosen for your mesh-refinement study converges with increasing polynomial order. This process can sometimes
be iterative; as one learns more about the flow and its physics, the need for additional h-refinement may become apparent after failing to achieve mesh convergence with p-refinement alone.

The use of higher-order meshes also allows for lower quality elements than typically seen in other finite-element solvers. Depending on the physics, an average aspect ratio of 20 or more can be tolerated for your mesh. Care must be taken to not create overly
small elements typical of low-order finite-element codes, since NekRS will allocate ``N+1`` Gauss-Lobatto Legendre (GLL) nodes per element, where ``N`` is the polynomial order defined in the ``.par`` file. This can easily lead to excessive degrees of freedom.
While the polynomial order can be reduced for such meshes (for example, due to small features that impose element size limitations on the entire mesh), NekRS performs best when operating at higher polynomial orders (``N``>=5).
