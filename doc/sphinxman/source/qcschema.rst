QCSchema in Psi4
================

Psi4 can read and write `QCSchema <https://molssi.github.io/QCElemental/dev/models.html>`_ inputs/results,
both through the command line and through the Python API.
This is useful when you want a stable, machine-friendly interface for single-point
energies, gradients, and related data. It is not appropriate for geometry optimizations,
multi-stage computations, or accessing detailed data like integrals.

Psi4 supports both QCSchema v1 (:py:func:`qcelemental.models.v1.AtomicInput` / :py:func:`qcelemental.models.v1.AtomicResult`) and
QCSchema v2 (:py:func:`qcelemental.models.v2.AtomicInput` / :py:func:`qcelemental.models.v2.AtomicResult`).


Command-Line Usage
------------------

Use ``--qcschema`` to run an input file as QCSchema:

.. code-block:: bash

   psi4 --qcschema input.json

If ``-o`` is omitted, Psi4 writes the result back to the same filename. To keep
input and output separate, provide ``-o``:

.. code-block:: bash

   psi4 --qcschema input.json -o output.json

You can also request a specific output schema version:

.. code-block:: bash

   psi4 --qcschema input.json --return-version 2 -o output.json


Minimal v2 Input (CLI or API)
-----------------------------

.. code-block:: json

   {
     "schema_name": "qcschema_atomic_input",
     "schema_version": 2,
     "molecule": {
       "symbols": ["O", "H", "H"],
       "geometry": [
         0.0, 0.0, -0.124,
         0.0, -1.432, 0.986,
         0.0,  1.432, 0.986
       ]
     },
     "specification": {
       "driver": "energy",
       "model": {"method": "scf", "basis": "cc-pVDZ"},
        "keywords": {"scf_type": "df"},
        "protocols": {"stdout": true}
      }
    }


Python API: ``run_qcschema``
----------------------------

Programmatic entry point is ``psi4.schema_wrapper.run_qcschema`` (also available as ``psi4.run_qcschema``).

.. autofunction:: psi4.schema_wrapper.run_qcschema

Model return (default):

.. code-block:: python

   import pprint
   import psi4

   inp = {
       "schema_name": "qcschema_atomic_input",
       "schema_version": 2,
       "molecule": {
           "symbols": ["He", "He"],
           "geometry": [0.0, 0.0, -2.5, 0.0, 0.0, 2.5]
       },
       "specification": {
           "driver": "energy",
           "model": {"method": "scf", "basis": "cc-pVDZ"},
           "keywords": {"scf_type": "df"}
       }
   }

   ret = psi4.run_qcschema(inp)
   print(ret.success)
   print(ret.return_result)
   pprint.pprint(ret.model_dump(), width=200)

Dictionary return:

.. code-block:: python

   ret = psi4.schema_wrapper.run_qcschema(inp, return_dict=True)
   print(ret["success"])
   print(ret["return_result"])
   print(ret["extras"]["qcvars"]["SCF TOTAL ENERGY"])


Serialization Notes
-------------------

- For v1 model objects, use ``.dict()`` for Python dictionary or ``.json()`` for JSON.
- For v2 model objects, use ``.model_dump()`` or ``.model_dump_json()``, respectively.
- If you request ``return_dict=True``, you already have plain Python containers
  and can use ``json.dump(...)`` directly. Note that with Python 3.14 or later, you
  can only use QCSchema v1 with dictionaries, not Pydantic models, due to Pydantic
  restrictions. See the ``run_qcschema`` docstring for details.


Running Through QCEngine
------------------------

QCEngine can run Psi4 as a backend using the same QCSchema payload.

.. code-block:: python

   import qcengine as qcng

   inp = {
       "schema_name": "qcschema_atomic_input",
       "schema_version": 2,
       "molecule": {
           "symbols": ["O", "H", "H"],
           "geometry": [0.0, 0.0, -0.124, 0.0, -1.432, 0.986, 0.0, 1.432, 0.986],
       },
       "specification": {
           "driver": "energy",
           "model": {"method": "scf", "basis": "cc-pVDZ"},
           "keywords": {"scf_type": "df"},
       },
   }

   ret = qcng.compute(inp, "psi4", return_dict=True)
   print(ret["success"])
   print(ret["return_result"])

You may also pass ``return_version=1`` or ``return_version=2`` to control
the output schema version.


Where to Find More Examples
---------------------------

- API-style integration tests:
  :source:`tests/json/schema-1-energy/input.py`,
  :source:`tests/json/schema-2-energy/input.py`,
  :source:`tests/json/schema-2-gradient/input.py`,
  :source:`tests/json/schema-2-properties/input.py`,
  :source:`tests/json/schema-2-throws/input.py`.
- API and CLI coverage in :source:`tests/pytests/test_qcschema.py`, especially
  ``test_qcschema_energy``, ``test_qcschema_gradient``, ``test_qcschema_cli``,
  and ``test_qcschema_wavefunction_basis``.
- Additional addon schema tests in :source:`tests/pytests/test_addons_qcschema.py`.
