.. _boundary_conditions:

Boundary conditions
===================

Boundary conditions for each mesh boundary should normally be set in the :ref:`parameter_file`, using the ``boundaryTypeMap`` parameter.
This is used within the ``VELOCITY``, ``TEMPERATURE`` or ``SCALARXX`` sections to set the boundary conditions of the respective solvers for the case.

Available Types
---------------

All available boundary conditions are summarised in the command line help function for nekRS (``nekrs --help par``) and below.

.. literalinclude:: ../../parHelp.txt
   :language: none
   :lines: 1-3, 149-180

.. tip::

   To setup some cases (e.g., cases that use third party meshes wherein boundaries are identified by a unique integer), you may need to use the ``boundaryIDMap`` parameter of the ``MESH`` section to apply the ``boundaryTypeMap`` options to the correct boundary IDs of the mesh (By default NekRS assumes that the boundaryIDs start at 1).
   Assuming an inlet-outlet pipe flow, below is an example where the three boundary conditions are applied to the boundary IDs 389 (``codeFixedValue`` - inlet), 231 (``zeroGradient`` - outlet), 4 (``zeroValue`` - walls).

   .. code-block::

      [VELOCITY]
      boundaryTypeMap = codedFixedValue, zeroGradient, zeroValue

      [MESH]
      boundaryIDMap = 389, 231, 4

.. note::

  For conjugate heat transfer cases, it is necessary to also specify the velocity boundary map, ``boundaryIDMapv``, separately to identify the interface conformal walls between the fluid and solid mesh.
  Following the above example, assume that the pipe is enclosed in a solid annular insulated jacket where the boundary ID 4 corresponds to the conformal wall between the solid and fluid mesh, while remaining insulated walls of the solid mesh are 10, 20 and 30.

    .. code-block::

       [MESH]
       boundaryIDMapV = codedFixedValue, zeroGradient, zeroValue
       boundaryIDMap = zeroGradient, zeroGradient, zeroGradient, codedFixedValue, zeroGradient

       [VELOCITY]
       boundaryTypeMap = 389, 231 ,4

       [TEMPERATURE]
       boundaryTypeMap = 10, 20, 30, 389, 231

User Defined Value/Gradient
"""""""""""""""""""""""""""

If a boundary condition requires a value setting rather than just a type, a suitable function will need to be provided within the ``.oudf`` file (or :ref:`okl block <okl_block>` section of ``.udf`` file).
The name of the function used should be in the form of ``codedFixed`` + ``Value/Gradient`` + ``Velocity/Scalar/Mesh``, E.G. ``codedFixedValueVelocity``, ``codedFixedGradientScalar`` and ``codedFixedValueMesh``.

All of these functions are passed the ``bcData`` struct which has the following parameters available:

+-----------------------------------------------------+------------------------------+---------------------------------------+
| Name                                                | Type                         | Description                           |
+=====================================================+==============================+=======================================+
| ``idM``                                             | ``int``                      | mesh index of boundary GLL            |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``fieldOffset``                                     | ``int``                      | array offset                          |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``id``                                              | ``int``                      | boundaryID                            |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``time``                                            | ``double``                   | current simulation time               |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``x``, ``y``, ``z``                                 | ``dfloat``                   | location coordinates of boundary GLL  |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``nx``, ``ny``, ``nz``                              | ``dfloat``                   | local normal vector components        |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``t1x``, ``t1y``, ``t1z``                           | ``dfloat``                   | local tangent vector components       |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``t2x``, ``t2y``, ``t2z``                           | ``dfloat``                   | local bi-tangent vector components    |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``tr1``                                             | ``dfloat``                   | local traction along tangent          |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``tr2``                                             | ``dfloat``                   | local traction along bi-tangent       |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``u``, ``v``, ``w``                                 | ``dfloat``                   | local velocity components             |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``p``                                               | ``dfloat``                   | local pressure                        |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``uinterp``, ``vinterp``, ``winterp``               | ``dfloat``                   | interpolated velocity components      |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``scalarId``                                        | ``int``                      | scalar index                          |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``s``                                               | ``dfloat``                   | scalar value (Dirichlet)              |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``flux``                                            | ``dfloat``                   | scalar flux (Neumann)                 |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``sinterp``                                         | ``dfloat``                   | interpolated scalar value             |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``meshu``, ``meshv``, ``meshw``                     | ``dfloat``                   | mesh velocity components              |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``trans``                                           | ``dfloat``                   | :math:`\rho` or :math:`\rho c_p`      |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``diff``                                            | ``dfloat``                   | diffusion coefficient                 |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``usrwrk``                                          | ``@globalPtr const dfloat*`` | pointer to user array                 |
+-----------------------------------------------------+------------------------------+---------------------------------------+


.. note::

    Some example cases are noted here for different user defined boundary condition types:

    +-----------------------------------------------------+------------------------------------------------------------------------------------+
    | BC Type                                             | Example Case                                                                       |                
    +=====================================================+====================================================================================+
    | ``codedFixedValueVelocity``                         | - `gabls <https://github.com/Nek5000/nekRS/tree/next/examples/gabls1>`_            |
    |                                                     | - `lowMach <https://github.com/Nek5000/nekRS/tree/next/examples/lowMach>`_         | 
    |                                                     | - `hemi <https://github.com/Nek5000/nekRS/tree/next/examples/hemi>`_               | 
    |                                                     | - `conj_ht <https://github.com/Nek5000/nekRS/tree/next/examples/conj_ht>`_         | 
    |                                                     | - `mv_cyl <https://github.com/Nek5000/nekRS/tree/next/examples/mv_cyl>`_           | 
    +-----------------------------------------------------+------------------------------------------------------------------------------------+
    | ``codedFixedValueVelocity`` (interpolated)          | `eddyNekNek <https://github.com/Nek5000/nekRS/tree/next/examples/eddyNekNek>`_     |
    +-----------------------------------------------------+------------------------------------------------------------------------------------+
    | ``codedFixedGradientVelocity``                      | `gabls <https://github.com/Nek5000/nekRS/tree/next/examples/gabls1>`_              |
    +-----------------------------------------------------+------------------------------------------------------------------------------------+
    | ``codedFixedGradientScalar``                        | `gabls <https://github.com/Nek5000/nekRS/tree/next/examples/gabls1>`_              |
    +-----------------------------------------------------+------------------------------------------------------------------------------------+
    | ``codedFixedValueScalar``                           | - `lowMach <https://github.com/Nek5000/nekRS/tree/next/examples/lowMach>`_         |
    |                                                     | - `conj_ht <https://github.com/Nek5000/nekRS/tree/next/examples/conj_ht>`_         | 
    |                                                     | - `ktauChannel <https://github.com/Nek5000/nekRS/tree/next/examples/ktauChannel>`_ | 
    +-----------------------------------------------------+------------------------------------------------------------------------------------+
    | ``codedFixedValueScalar`` (interpolated)            | `eddyNekNek <https://github.com/Nek5000/nekRS/tree/next/examples/eddyNekNek>`_     |
    +-----------------------------------------------------+------------------------------------------------------------------------------------+
    | ``codedFixedValuePressure``                         | `hemi <https://github.com/Nek5000/nekRS/tree/next/examples/hemi>`_                 |
    +-----------------------------------------------------+------------------------------------------------------------------------------------+

    

Internal / Periodic
"""""""""""""""""""
.. _periodic_boundary:

``None`` is used when a internal boundary condition is required or a periodic boundary condition has been set as part of the mesh and it does not need to be considered as part of the standard processing of boundary conditions.
Perodicity is linked to the mesh connectivity and is handled by the meshing tool.


Turbulent Inflow (for LES/DNS)
----------------------

Velocity Recycling
""""""""""""""""""

NekRS offers an in-built technique to impose fully developed turbulent inflow.  
It comprises recycling the velocity from an offset surface in the domain interior to the inlet boundary resulting in a quasi-periodic domain. 
The required routines are available through the ``velRecycling.hpp`` header file. 
Sample ``.udf`` code to use velocity recycling is as follows,

.. code-block::

  #include "velRecycling.hpp"

  #ifdef __okl__
   void codedFixedValueVelocity(bcData * bc)
   {
      bc->u = bc->usrwrk[bc->idM + 0 * bc->fieldOffset];
      bc->v = bc->usrwrk[bc->idM + 1 * bc->fieldOffset];
      bc->w = bc->usrwrk[bc->idM + 2 * bc->fieldOffset];
   }
  #endif

  void UDF_Setup()
  {
    nrs->o_usrwrk.resize(nrs->NVfields * nrs->fieldOffset);

    const dfloat uBulk = 1.0;                              //mean bulk velocity
    const int bID = 1;                                     //boundary ID of inlet
    const dfloat Dh = 1.0;                                 //Hydraulic diameter
    condt dfloat xOffset = 0.0;                            //x-offset of recycling surface
    const dfloat yOffset = 0.0;                            //y-offset of recycling surface
    const dfloat zOffset = 10.0 * Dh;                      //z-offset of recycling surface
    velRecycling::setup(nrs->o_usrwrk, xOffset, yOffset, zOffset, bID, uBulk);
  }
  void UDF_ExecuteStep(double time, int tstep)
  {
    velRecycling::copy();
  }

The above example assumes that the inlet surface, with boundary ID (``bID``) equal to 1, is perpendicular to the *x-y* plane.
The recycling surface is located at an offset distance equal to :math:`10 Dh` from the inlet surface, where :math:`Dh` is the hydraulic diameter of the inlet cross-section. 
``uBulk`` is the desired mean velocity at the inlet boundary. 

.. note::

  The recommended recycling length is :math:`8-10 Dh` [Lund1998]_ to allow the resolution of largest eddies. The user must ensure that the boundary layer is not acted upon by significant pressure gradients or geometrical changes between the inlet and recycle plane.

The interpolated velocities from the recycling surface are stored in the ``nrs->o_usrwrk`` array which is accessible in the ``codedFixedValueVelocity`` kernel.
Therefore, ``nrs->o_usrwrk`` vector must be resized to accomodate the velocity data to ``nrs->NVfields * nrs->fieldOffset`` length.
The ``velRecycling::setup`` call setups up the necessary interpolation procedures to map the velocity from domain interior to the inlet surface.
Note the parameters passed as arguments to the setup routine.
The ``velRecycling::copy`` routine is called from the ``UDF_ExecuteStep`` routine at every time step.
It, internally, populates the ``nrs->o_usrwrk`` array with the scaled instantaneous velocity at the recycle plane based on specified ``ubulk``.
Finally, the captured velocity is specified as inlet Dirichlet boundary condition in ``codedFixedValueVelocity`` in the :ref:`okl block <okl_block>` section, as shown above.
Note the proper indexing and offset specified using ``bc->idM`` and ``bc->fieldOffset`` for the three velocity components.
