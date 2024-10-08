.. _postprocessing:

Postprocessing
==============

Once a case has been setup correctly so it can be run without errors, you may want 
to modify the postprocessing to of the simulation output. By default, NekRS will
output a basic set of data according to frequency set in the ``writeInterval`` of
the :ref:`parameter_file` which can subsequently be viewed through a visualization
tool such as Paraview or Visit. However, additional data or derived values can
be extracted by setting up User Defined outputs using the ``UDF_ExecuteStep``
function of the ``.udf`` file.

.. _checkpointing_visualisation:

Checkpointing & Visualization
-----------------------------

Standard NekRS field data output files have the form ``<case>0.f<n>``, where ``<case>`` is the
case name and ``<n>`` is a five-digit number indicating the number of the output
file (each output file represents a single time step that is output according to
the settings for ``writeControl`` and ``writeInterval`` in the ``.par`` file).
These output files are in a custom binary format that requires a ``<case>.nek5000``
file to be viewable in Paraview (or similar tools). This should be automatically
generated by NekRS, but can also be manually created using the ``nrsvis`` script.

.. _compute_derived:

Compute Derived Quantity
------------------------

Additional control of the simulation to compute additional/derived quantities 
or output custom fields can be achieved by utilising the ``UDF_ExecuteStep`` 
function of the ``.udf`` file. Here we demonstrate how this can be used to 
compute a derived quantity and output custom fields.

Qcriterion in turbPipe example

.. _custom_checkpoint:

Adding Custom Checkpoint Fields
-------------------------------

.. code-block:: cpp
    
    nrs->addUserCheckpointField("scalar01", std::vector<deviceMemory<dfloat>>{o_nuAVM});

Adding Custom Output File
-------------------------

.. code-block:: cpp

    iofld = iofldFactory::create();
    iofld->open(mesh, iofld::mode::write, "qcrit");
    iofld->writeAttribute("uniform", "true");
    iofld->writeAttribute("polynomialOrder", std::to_string(mesh->N + 2));
    iofld->process();

.. _turbulence_stats:

Run time averaging
------------------
