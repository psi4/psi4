#! test QCSchema for response properties

import psi4
import numpy as np
import json
import os


# Generate JSON input
json_data = {
    "schema_name": "qcschema_input",
    "schema_version": 1,
    "molecule": {
        "geometry": [
            0.0, 0.0, -0.1294769411935893, 0.0, -1.494187339479985, 1.0274465079245698, 0.0, 1.494187339479985,
            1.0274465079245698
        ],
        "symbols": ["O", "H", "H"],
        "fix_com":
        True,
        "fix_orientation":
        True
    },
    "driver": "properties",
    "model": {
        "method": "CC2",
        "basis": "6-31G",
    },
    "keywords": {
        "scf_type": "df",
        "mp2_type": "df",
        "e_convergence": 9,
        "omega": [355, 439, 'nm'],
        "gauge": "velocity",
        "function_kwargs": {
            "properties": ["dipole", "polarizability", "rotation", "roa_tensor"]
        }
    }
}

# Load expected output (dipole & quadrupole in au)
expected_response = {
  "CC2 DIPOLE POLARIZABILITY @ 355NM": 5.175390333179149,
  "CC2 DIPOLE POLARIZABILITY @ 439NM": 5.063282247200858,
  #"CC2 DIPOLE X": 0.0,
  #"CC2 DIPOLE Y": 0.0,
  #"CC2 DIPOLE Z": 2.5171476216528084,
  #"CC2 QUADRUPOLE XX": -7.393278028527795,
  #"CC2 QUADRUPOLE XY": -1.8387569175060344e-16,
  #"CC2 QUADRUPOLE XZ": 0.0,
  #"CC2 QUADRUPOLE YY": -4.3142280005261755,
  #"CC2 QUADRUPOLE YZ": 0.0,
  #"CC2 QUADRUPOLE ZZ": -5.9972667929848065,
  "CC2 SPECIFIC ROTATION (MVG) @ 355NM": -0.0,
  "CC2 SPECIFIC ROTATION (MVG) @ 439NM": -0.0,
  "CC2 SPECIFIC ROTATION (VEL) @ 355NM": -0.0,
  "CC2 SPECIFIC ROTATION (VEL) @ 439NM": -0.0,
  "CC2 DIPOLE": [
    0,
    0,
    0.99032208
  ],
  "CC2 QUADRUPOLE": [
    [
      -5.49672082,
      -1.36707065e-16,
      0.0
    ],
    [
      -1.36707065e-16,
      -3.20752267,
      0.0
    ],
    [
      0.0,
      0.0,
      -4.45882072
    ]
  ],
  "CC2 DIPOLE POLARIZABILITY TENSOR @ 355NM": [
    [
      1.53032264,
      0.0,
      0.0
    ],
    [
      0.0,
      8.37432166,
      0.0
    ],
    [
      0.0,
      0.0,
      5.62152674
    ]
  ],
  "CC2 DIPOLE POLARIZABILITY TENSOR @ 439NM": [
    [
      1.49010149,
      0.0,
      0.0
    ],
    [
      0.0,
      8.21869196,
      0.0
    ],
    [
      0.0,
      0.0,
      5.48105331
    ]
  ],
  "CC2 OPTICAL ROTATION TENSOR (MVG) @ 355NM": [
    [
      0.0,
      -0.02923246,
      0.0
    ],
    [
      -0.01267525,
      0.0,
      0.0
    ],
    [
      0.0,
      0.0,
      0.0
    ]
  ],
  "CC2 OPTICAL ROTATION TENSOR (MVG) @ 439NM": [
    [
      0.0,
      -0.01758964,
      0.0
    ],
    [
      -0.00805598,
      0.0,
      0.0
    ],
    [
      0.0,
      0.0,
      0.0
    ]
  ],
  "CC2 OPTICAL ROTATION TENSOR (VEL) @ 0NM": [
    [
      0.0,
      0.0124954,
      0.0
    ],
    [
      -0.25085903,
      0.0,
      0.0
    ],
    [
      0.0,
      0.0,
      0.0
    ]
  ],
  "CC2 OPTICAL ROTATION TENSOR (VEL) @ 355NM": [
    [
      0.0,
      -0.01673706,
      0.0
    ],
    [
      -0.26353427,
      0.0,
      0.0
    ],
    [
      0.0,
      0.0,
      0.0
    ]
  ],
  "CC2 OPTICAL ROTATION TENSOR (VEL) @ 439NM": [
    [
      0.0,
      -0.00509424,
      0.0
    ],
    [
      -0.25891501,
      0.0,
      0.0
    ],
    [
      0.0,
      0.0,
      0.0
    ]
  ],
  "CC2 QUADRUPOLE POLARIZABILITY TENSOR @ 355NM": [[
    [
      0.0,
      0.0,
      0.38619772
    ],
    [
      0.0,
      0.0,
      0.0
    ],
    [
      0.38619772,
      0.0,
      0.0
    ]
  ], [
    [
      0.0,
      0.0,
      0.0
    ],
    [
      0.0,
      0.0,
      10.97015701
    ],
    [
      0.0,
      10.97015701,
      0.0
    ]
  ], [
    [
      -5.65001748,
      0.0,
      0.0
    ],
    [
      0.0,
      5.62557634,
      0.0
    ],
    [
      0.0,
      0.0,
      0.02444115
    ]
  ]],
  "CC2 QUADRUPOLE POLARIZABILITY TENSOR @ 439NM": [[
    [
      0.0,
      0.0,
      0.37061096
    ],
    [
      0.0,
      0.0,
      0.0
    ],
    [
      0.37061096,
      0.0,
      0.0
    ]
  ], [
    [
      0.0,
      0.0,
      0.0
    ],
    [
      0.0,
      0.0,
      10.74076995
    ],
    [
      0.0,
      10.74076995,
      0.0
    ]
  ], [
    [
      -5.51540341,
      0.0,
      0.0
    ],
    [
      0.0,
      5.49921728,
      0.0
    ],
    [
      0.0,
      0.0,
      0.01618613
    ]
  ]]
}

expected_response = {k: (np.asarray(v) if isinstance(v, list) else v) for k, v in expected_response.items()}

json_ret = psi4.schema_wrapper.run_qcschema(json_data).dict()

