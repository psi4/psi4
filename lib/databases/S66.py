"""
| Database (Hobza) of interaction energies for bimolecular complexes.
| Geometries and reference energies from Rezac et al. JCTC 7 2427 (2011).

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'``

- **subset**

  - ``'small'``
  - ``'large'``
  - ``'HB'`` hydrogen-bonded systems
  - ``'MX'`` mixed-influence systems
  - ``'DD'`` dispersion-dominated systems

"""
import re
import input

# <<< S66 Database Module >>>
dbse = 'S66'

# <<< Database Members >>>
HRXN = range(1, 67)
HRXN_SM = [1, 12, 59]
HRXN_LG = [26, 34]
HB = range(1, 24)
MX = range(47, 67)
DD = range(24, 47)

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}  # order of active reagents for non-supermolecular calculations
for rxn in HRXN:

    RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                      '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                      '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                      '%s-%s-monoA-unCP' % (dbse, rxn) : -1,
                                      '%s-%s-monoB-unCP' % (dbse, rxn) : -1 }

    ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

    ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                      '%s-%s-monoA-CP'   % (dbse, rxn),
                                      '%s-%s-monoB-CP'   % (dbse, rxn) ]

    ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                      '%s-%s-monoA-unCP' % (dbse, rxn),
                                      '%s-%s-monoB-unCP' % (dbse, rxn) ]

# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'            % (dbse, '1'                     )] =   -4.918
BIND['%s-%s'            % (dbse, '2'                     )] =   -5.592
BIND['%s-%s'            % (dbse, '3'                     )] =   -6.908
BIND['%s-%s'            % (dbse, '4'                     )] =   -8.103
BIND['%s-%s'            % (dbse, '5'                     )] =   -5.757
BIND['%s-%s'            % (dbse, '6'                     )] =   -7.554
BIND['%s-%s'            % (dbse, '7'                     )] =   -8.230
BIND['%s-%s'            % (dbse, '8'                     )] =   -5.009
BIND['%s-%s'            % (dbse, '9'                     )] =   -3.059
BIND['%s-%s'            % (dbse, '10'                    )] =   -4.160
BIND['%s-%s'            % (dbse, '11'                    )] =   -5.419
BIND['%s-%s'            % (dbse, '12'                    )] =   -7.266
BIND['%s-%s'            % (dbse, '13'                    )] =   -6.187
BIND['%s-%s'            % (dbse, '14'                    )] =   -7.454
BIND['%s-%s'            % (dbse, '15'                    )] =   -8.630
BIND['%s-%s'            % (dbse, '16'                    )] =   -5.124
BIND['%s-%s'            % (dbse, '17'                    )] =  -17.182
BIND['%s-%s'            % (dbse, '18'                    )] =   -6.857
BIND['%s-%s'            % (dbse, '19'                    )] =   -7.410
BIND['%s-%s'            % (dbse, '20'                    )] =  -19.093
BIND['%s-%s'            % (dbse, '21'                    )] =  -16.265
BIND['%s-%s'            % (dbse, '22'                    )] =  -19.491
BIND['%s-%s'            % (dbse, '23'                    )] =  -19.189
BIND['%s-%s'            % (dbse, '24'                    )] =   -2.822
BIND['%s-%s'            % (dbse, '25'                    )] =   -3.895
BIND['%s-%s'            % (dbse, '26'                    )] =   -9.829
BIND['%s-%s'            % (dbse, '27'                    )] =   -3.439
BIND['%s-%s'            % (dbse, '28'                    )] =   -5.713
BIND['%s-%s'            % (dbse, '29'                    )] =   -6.819
BIND['%s-%s'            % (dbse, '30'                    )] =   -1.432
BIND['%s-%s'            % (dbse, '31'                    )] =   -3.380
BIND['%s-%s'            % (dbse, '32'                    )] =   -3.738
BIND['%s-%s'            % (dbse, '33'                    )] =   -1.872
BIND['%s-%s'            % (dbse, '34'                    )] =   -3.776
BIND['%s-%s'            % (dbse, '35'                    )] =   -2.613
BIND['%s-%s'            % (dbse, '36'                    )] =   -1.777
BIND['%s-%s'            % (dbse, '37'                    )] =   -2.404
BIND['%s-%s'            % (dbse, '38'                    )] =   -2.997
BIND['%s-%s'            % (dbse, '39'                    )] =   -3.575
BIND['%s-%s'            % (dbse, '40'                    )] =   -2.895
BIND['%s-%s'            % (dbse, '41'                    )] =   -4.848
BIND['%s-%s'            % (dbse, '42'                    )] =   -4.138
BIND['%s-%s'            % (dbse, '43'                    )] =   -3.712
BIND['%s-%s'            % (dbse, '44'                    )] =   -2.005
BIND['%s-%s'            % (dbse, '45'                    )] =   -1.748
BIND['%s-%s'            % (dbse, '46'                    )] =   -4.264
BIND['%s-%s'            % (dbse, '47'                    )] =   -2.876
BIND['%s-%s'            % (dbse, '48'                    )] =   -3.535
BIND['%s-%s'            % (dbse, '49'                    )] =   -3.331
BIND['%s-%s'            % (dbse, '50'                    )] =   -2.867
BIND['%s-%s'            % (dbse, '51'                    )] =   -1.524
BIND['%s-%s'            % (dbse, '52'                    )] =   -4.707
BIND['%s-%s'            % (dbse, '53'                    )] =   -4.361
BIND['%s-%s'            % (dbse, '54'                    )] =   -3.277
BIND['%s-%s'            % (dbse, '55'                    )] =   -4.188
BIND['%s-%s'            % (dbse, '56'                    )] =   -3.231
BIND['%s-%s'            % (dbse, '57'                    )] =   -5.282
BIND['%s-%s'            % (dbse, '58'                    )] =   -4.146
BIND['%s-%s'            % (dbse, '59'                    )] =   -2.850
BIND['%s-%s'            % (dbse, '60'                    )] =   -4.868
BIND['%s-%s'            % (dbse, '61'                    )] =   -2.912
BIND['%s-%s'            % (dbse, '62'                    )] =   -3.534
BIND['%s-%s'            % (dbse, '63'                    )] =   -3.801
BIND['%s-%s'            % (dbse, '64'                    )] =   -2.999
BIND['%s-%s'            % (dbse, '65'                    )] =   -3.991
BIND['%s-%s'            % (dbse, '66'                    )] =   -3.968

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, '1'                     )] = """Water Dimer """
TAGL['%s-%s-dimer'      % (dbse, '1'                     )] = """Dimer from Water Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '1'                     )] = """Monomer A from Water Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '1'                     )] = """Monomer B from Water Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '1'                     )] = """Monomer A from Water Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '1'                     )] = """Monomer B from Water Dimer """
TAGL['%s-%s'            % (dbse, '2'                     )] = """Water-Methanol """
TAGL['%s-%s-dimer'      % (dbse, '2'                     )] = """Dimer from Water-Methanol """
TAGL['%s-%s-monoA-CP'   % (dbse, '2'                     )] = """Monomer A from Water-Methanol """
TAGL['%s-%s-monoB-CP'   % (dbse, '2'                     )] = """Monomer B from Water-Methanol """
TAGL['%s-%s-monoA-unCP' % (dbse, '2'                     )] = """Monomer A from Water-Methanol """
TAGL['%s-%s-monoB-unCP' % (dbse, '2'                     )] = """Monomer B from Water-Methanol """
TAGL['%s-%s'            % (dbse, '3'                     )] = """Water-Methylamine """
TAGL['%s-%s-dimer'      % (dbse, '3'                     )] = """Dimer from Water-Methylamine """
TAGL['%s-%s-monoA-CP'   % (dbse, '3'                     )] = """Monomer A from Water-Methylamine """
TAGL['%s-%s-monoB-CP'   % (dbse, '3'                     )] = """Monomer B from Water-Methylamine """
TAGL['%s-%s-monoA-unCP' % (dbse, '3'                     )] = """Monomer A from Water-Methylamine """
TAGL['%s-%s-monoB-unCP' % (dbse, '3'                     )] = """Monomer B from Water-Methylamine """
TAGL['%s-%s'            % (dbse, '4'                     )] = """Water-N-methylacetamide """
TAGL['%s-%s-dimer'      % (dbse, '4'                     )] = """Dimer from Water-N-methylacetamide """
TAGL['%s-%s-monoA-CP'   % (dbse, '4'                     )] = """Monomer A from Water-N-methylacetamide """
TAGL['%s-%s-monoB-CP'   % (dbse, '4'                     )] = """Monomer B from Water-N-methylacetamide """
TAGL['%s-%s-monoA-unCP' % (dbse, '4'                     )] = """Monomer A from Water-N-methylacetamide """
TAGL['%s-%s-monoB-unCP' % (dbse, '4'                     )] = """Monomer B from Water-N-methylacetamide """
TAGL['%s-%s'            % (dbse, '5'                     )] = """Methanol Dimer """
TAGL['%s-%s-dimer'      % (dbse, '5'                     )] = """Dimer from Methanol Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '5'                     )] = """Monomer A from Methanol Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '5'                     )] = """Monomer B from Methanol Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '5'                     )] = """Monomer A from Methanol Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '5'                     )] = """Monomer B from Methanol Dimer """
TAGL['%s-%s'            % (dbse, '6'                     )] = """Methanol-Methylamine """
TAGL['%s-%s-dimer'      % (dbse, '6'                     )] = """Dimer from Methanol-Methylamine """
TAGL['%s-%s-monoA-CP'   % (dbse, '6'                     )] = """Monomer A from Methanol-Methylamine """
TAGL['%s-%s-monoB-CP'   % (dbse, '6'                     )] = """Monomer B from Methanol-Methylamine """
TAGL['%s-%s-monoA-unCP' % (dbse, '6'                     )] = """Monomer A from Methanol-Methylamine """
TAGL['%s-%s-monoB-unCP' % (dbse, '6'                     )] = """Monomer B from Methanol-Methylamine """
TAGL['%s-%s'            % (dbse, '7'                     )] = """Methanol-N-methylacetamide """
TAGL['%s-%s-dimer'      % (dbse, '7'                     )] = """Dimer from Methanol-N-methylacetamide """
TAGL['%s-%s-monoA-CP'   % (dbse, '7'                     )] = """Monomer A from Methanol-N-methylacetamide """
TAGL['%s-%s-monoB-CP'   % (dbse, '7'                     )] = """Monomer B from Methanol-N-methylacetamide """
TAGL['%s-%s-monoA-unCP' % (dbse, '7'                     )] = """Monomer A from Methanol-N-methylacetamide """
TAGL['%s-%s-monoB-unCP' % (dbse, '7'                     )] = """Monomer B from Methanol-N-methylacetamide """
TAGL['%s-%s'            % (dbse, '8'                     )] = """Methanol-Water """
TAGL['%s-%s-dimer'      % (dbse, '8'                     )] = """Dimer from Methanol-Water """
TAGL['%s-%s-monoA-CP'   % (dbse, '8'                     )] = """Monomer A from Methanol-Water """
TAGL['%s-%s-monoB-CP'   % (dbse, '8'                     )] = """Monomer B from Methanol-Water """
TAGL['%s-%s-monoA-unCP' % (dbse, '8'                     )] = """Monomer A from Methanol-Water """
TAGL['%s-%s-monoB-unCP' % (dbse, '8'                     )] = """Monomer B from Methanol-Water """
TAGL['%s-%s'            % (dbse, '9'                     )] = """Methylamine-Methanol """
TAGL['%s-%s-dimer'      % (dbse, '9'                     )] = """Dimer from Methylamine-Methanol """
TAGL['%s-%s-monoA-CP'   % (dbse, '9'                     )] = """Monomer A from Methylamine-Methanol """
TAGL['%s-%s-monoB-CP'   % (dbse, '9'                     )] = """Monomer B from Methylamine-Methanol """
TAGL['%s-%s-monoA-unCP' % (dbse, '9'                     )] = """Monomer A from Methylamine-Methanol """
TAGL['%s-%s-monoB-unCP' % (dbse, '9'                     )] = """Monomer B from Methylamine-Methanol """
TAGL['%s-%s'            % (dbse, '10'                    )] = """Methylamine Dimer """
TAGL['%s-%s-dimer'      % (dbse, '10'                    )] = """Dimer from Methylamine Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '10'                    )] = """Monomer A from Methylamine Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '10'                    )] = """Monomer B from Methylamine Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '10'                    )] = """Monomer A from Methylamine Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '10'                    )] = """Monomer B from Methylamine Dimer """
TAGL['%s-%s'            % (dbse, '11'                    )] = """Methylamine-N-methylacetamide """
TAGL['%s-%s-dimer'      % (dbse, '11'                    )] = """Dimer from Methylamine-N-methylacetamide """
TAGL['%s-%s-monoA-CP'   % (dbse, '11'                    )] = """Monomer A from Methylamine-N-methylacetamide """
TAGL['%s-%s-monoB-CP'   % (dbse, '11'                    )] = """Monomer B from Methylamine-N-methylacetamide """
TAGL['%s-%s-monoA-unCP' % (dbse, '11'                    )] = """Monomer A from Methylamine-N-methylacetamide """
TAGL['%s-%s-monoB-unCP' % (dbse, '11'                    )] = """Monomer B from Methylamine-N-methylacetamide """
TAGL['%s-%s'            % (dbse, '12'                    )] = """Methylamine-Water """
TAGL['%s-%s-dimer'      % (dbse, '12'                    )] = """Dimer from Methylamine-Water """
TAGL['%s-%s-monoA-CP'   % (dbse, '12'                    )] = """Monomer A from Methylamine-Water """
TAGL['%s-%s-monoB-CP'   % (dbse, '12'                    )] = """Monomer B from Methylamine-Water """
TAGL['%s-%s-monoA-unCP' % (dbse, '12'                    )] = """Monomer A from Methylamine-Water """
TAGL['%s-%s-monoB-unCP' % (dbse, '12'                    )] = """Monomer B from Methylamine-Water """
TAGL['%s-%s'            % (dbse, '13'                    )] = """N-methylacetamide-Methanol """
TAGL['%s-%s-dimer'      % (dbse, '13'                    )] = """Dimer from N-methylacetamide-Methanol """
TAGL['%s-%s-monoA-CP'   % (dbse, '13'                    )] = """Monomer A from N-methylacetamide-Methanol """
TAGL['%s-%s-monoB-CP'   % (dbse, '13'                    )] = """Monomer B from N-methylacetamide-Methanol """
TAGL['%s-%s-monoA-unCP' % (dbse, '13'                    )] = """Monomer A from N-methylacetamide-Methanol """
TAGL['%s-%s-monoB-unCP' % (dbse, '13'                    )] = """Monomer B from N-methylacetamide-Methanol """
TAGL['%s-%s'            % (dbse, '14'                    )] = """N-methylacetamide-Methylamine """
TAGL['%s-%s-dimer'      % (dbse, '14'                    )] = """Dimer from N-methylacetamide-Methylamine """
TAGL['%s-%s-monoA-CP'   % (dbse, '14'                    )] = """Monomer A from N-methylacetamide-Methylamine """
TAGL['%s-%s-monoB-CP'   % (dbse, '14'                    )] = """Monomer B from N-methylacetamide-Methylamine """
TAGL['%s-%s-monoA-unCP' % (dbse, '14'                    )] = """Monomer A from N-methylacetamide-Methylamine """
TAGL['%s-%s-monoB-unCP' % (dbse, '14'                    )] = """Monomer B from N-methylacetamide-Methylamine """
TAGL['%s-%s'            % (dbse, '15'                    )] = """N-methylacetamide Dimer """
TAGL['%s-%s-dimer'      % (dbse, '15'                    )] = """Dimer from N-methylacetamide Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '15'                    )] = """Monomer A from N-methylacetamide Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '15'                    )] = """Monomer B from N-methylacetamide Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '15'                    )] = """Monomer A from N-methylacetamide Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '15'                    )] = """Monomer B from N-methylacetamide Dimer """
TAGL['%s-%s'            % (dbse, '16'                    )] = """N-methylacetamide-Water """
TAGL['%s-%s-dimer'      % (dbse, '16'                    )] = """Dimer from N-methylacetamide-Water """
TAGL['%s-%s-monoA-CP'   % (dbse, '16'                    )] = """Monomer A from N-methylacetamide-Water """
TAGL['%s-%s-monoB-CP'   % (dbse, '16'                    )] = """Monomer B from N-methylacetamide-Water """
TAGL['%s-%s-monoA-unCP' % (dbse, '16'                    )] = """Monomer A from N-methylacetamide-Water """
TAGL['%s-%s-monoB-unCP' % (dbse, '16'                    )] = """Monomer B from N-methylacetamide-Water """
TAGL['%s-%s'            % (dbse, '17'                    )] = """Uracil Dimer, HB """
TAGL['%s-%s-dimer'      % (dbse, '17'                    )] = """Dimer from Uracil Dimer, HB """
TAGL['%s-%s-monoA-CP'   % (dbse, '17'                    )] = """Monomer A from Uracil Dimer, HB """
TAGL['%s-%s-monoB-CP'   % (dbse, '17'                    )] = """Monomer B from Uracil Dimer, HB """
TAGL['%s-%s-monoA-unCP' % (dbse, '17'                    )] = """Monomer A from Uracil Dimer, HB """
TAGL['%s-%s-monoB-unCP' % (dbse, '17'                    )] = """Monomer B from Uracil Dimer, HB """
TAGL['%s-%s'            % (dbse, '18'                    )] = """Water-Pyridine """
TAGL['%s-%s-dimer'      % (dbse, '18'                    )] = """Dimer from Water-Pyridine """
TAGL['%s-%s-monoA-CP'   % (dbse, '18'                    )] = """Monomer A from Water-Pyridine """
TAGL['%s-%s-monoB-CP'   % (dbse, '18'                    )] = """Monomer B from Water-Pyridine """
TAGL['%s-%s-monoA-unCP' % (dbse, '18'                    )] = """Monomer A from Water-Pyridine """
TAGL['%s-%s-monoB-unCP' % (dbse, '18'                    )] = """Monomer B from Water-Pyridine """
TAGL['%s-%s'            % (dbse, '19'                    )] = """Methanol-Pyridine """
TAGL['%s-%s-dimer'      % (dbse, '19'                    )] = """Dimer from Methanol-Pyridine """
TAGL['%s-%s-monoA-CP'   % (dbse, '19'                    )] = """Monomer A from Methanol-Pyridine """
TAGL['%s-%s-monoB-CP'   % (dbse, '19'                    )] = """Monomer B from Methanol-Pyridine """
TAGL['%s-%s-monoA-unCP' % (dbse, '19'                    )] = """Monomer A from Methanol-Pyridine """
TAGL['%s-%s-monoB-unCP' % (dbse, '19'                    )] = """Monomer B from Methanol-Pyridine """
TAGL['%s-%s'            % (dbse, '20'                    )] = """Acetic Acid Dimer """
TAGL['%s-%s-dimer'      % (dbse, '20'                    )] = """Dimer from Acetic Acid Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '20'                    )] = """Monomer A from Acetic Acid Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '20'                    )] = """Monomer B from Acetic Acid Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '20'                    )] = """Monomer A from Acetic Acid Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '20'                    )] = """Monomer B from Acetic Acid Dimer """
TAGL['%s-%s'            % (dbse, '21'                    )] = """Acetamide Dimer """
TAGL['%s-%s-dimer'      % (dbse, '21'                    )] = """Dimer from Acetamide Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '21'                    )] = """Monomer A from Acetamide Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '21'                    )] = """Monomer B from Acetamide Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '21'                    )] = """Monomer A from Acetamide Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '21'                    )] = """Monomer B from Acetamide Dimer """
TAGL['%s-%s'            % (dbse, '22'                    )] = """Acetic Acid-Uracil """
TAGL['%s-%s-dimer'      % (dbse, '22'                    )] = """Dimer from Acetic Acid-Uracil """
TAGL['%s-%s-monoA-CP'   % (dbse, '22'                    )] = """Monomer A from Acetic Acid-Uracil """
TAGL['%s-%s-monoB-CP'   % (dbse, '22'                    )] = """Monomer B from Acetic Acid-Uracil """
TAGL['%s-%s-monoA-unCP' % (dbse, '22'                    )] = """Monomer A from Acetic Acid-Uracil """
TAGL['%s-%s-monoB-unCP' % (dbse, '22'                    )] = """Monomer B from Acetic Acid-Uracil """
TAGL['%s-%s'            % (dbse, '23'                    )] = """Acetamide-Uracil """
TAGL['%s-%s-dimer'      % (dbse, '23'                    )] = """Dimer from Acetamide-Uracil """
TAGL['%s-%s-monoA-CP'   % (dbse, '23'                    )] = """Monomer A from Acetamide-Uracil """
TAGL['%s-%s-monoB-CP'   % (dbse, '23'                    )] = """Monomer B from Acetamide-Uracil """
TAGL['%s-%s-monoA-unCP' % (dbse, '23'                    )] = """Monomer A from Acetamide-Uracil """
TAGL['%s-%s-monoB-unCP' % (dbse, '23'                    )] = """Monomer B from Acetamide-Uracil """
TAGL['%s-%s'            % (dbse, '24'                    )] = """Benzene Dimer, pi-pi """
TAGL['%s-%s-dimer'      % (dbse, '24'                    )] = """Dimer from Benzene Dimer, pi-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '24'                    )] = """Monomer A from Benzene Dimer, pi-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '24'                    )] = """Monomer B from Benzene Dimer, pi-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '24'                    )] = """Monomer A from Benzene Dimer, pi-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '24'                    )] = """Monomer B from Benzene Dimer, pi-pi """
TAGL['%s-%s'            % (dbse, '25'                    )] = """Pyridine Dimer, pi-pi """
TAGL['%s-%s-dimer'      % (dbse, '25'                    )] = """Dimer from Pyridine Dimer, pi-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '25'                    )] = """Monomer A from Pyridine Dimer, pi-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '25'                    )] = """Monomer B from Pyridine Dimer, pi-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '25'                    )] = """Monomer A from Pyridine Dimer, pi-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '25'                    )] = """Monomer B from Pyridine Dimer, pi-pi """
TAGL['%s-%s'            % (dbse, '26'                    )] = """Uracil Dimer, pi-pi """
TAGL['%s-%s-dimer'      % (dbse, '26'                    )] = """Dimer from Uracil Dimer, pi-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '26'                    )] = """Monomer A from Uracil Dimer, pi-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '26'                    )] = """Monomer B from Uracil Dimer, pi-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '26'                    )] = """Monomer A from Uracil Dimer, pi-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '26'                    )] = """Monomer B from Uracil Dimer, pi-pi """
TAGL['%s-%s'            % (dbse, '27'                    )] = """Benzene-Pyridine, pi-pi """
TAGL['%s-%s-dimer'      % (dbse, '27'                    )] = """Dimer from Benzene-Pyridine, pi-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '27'                    )] = """Monomer A from Benzene-Pyridine, pi-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '27'                    )] = """Monomer B from Benzene-Pyridine, pi-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '27'                    )] = """Monomer A from Benzene-Pyridine, pi-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '27'                    )] = """Monomer B from Benzene-Pyridine, pi-pi """
TAGL['%s-%s'            % (dbse, '28'                    )] = """Benzene-Uracil, pi-pi """
TAGL['%s-%s-dimer'      % (dbse, '28'                    )] = """Dimer from Benzene-Uracil, pi-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '28'                    )] = """Monomer A from Benzene-Uracil, pi-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '28'                    )] = """Monomer B from Benzene-Uracil, pi-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '28'                    )] = """Monomer A from Benzene-Uracil, pi-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '28'                    )] = """Monomer B from Benzene-Uracil, pi-pi """
TAGL['%s-%s'            % (dbse, '29'                    )] = """Pyridine-Uracil, pi-pi """
TAGL['%s-%s-dimer'      % (dbse, '29'                    )] = """Dimer from Pyridine-Uracil, pi-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '29'                    )] = """Monomer A from Pyridine-Uracil, pi-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '29'                    )] = """Monomer B from Pyridine-Uracil, pi-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '29'                    )] = """Monomer A from Pyridine-Uracil, pi-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '29'                    )] = """Monomer B from Pyridine-Uracil, pi-pi """
TAGL['%s-%s'            % (dbse, '30'                    )] = """Benzene-Ethene """
TAGL['%s-%s-dimer'      % (dbse, '30'                    )] = """Dimer from Benzene-Ethene """
TAGL['%s-%s-monoA-CP'   % (dbse, '30'                    )] = """Monomer A from Benzene-Ethene """
TAGL['%s-%s-monoB-CP'   % (dbse, '30'                    )] = """Monomer B from Benzene-Ethene """
TAGL['%s-%s-monoA-unCP' % (dbse, '30'                    )] = """Monomer A from Benzene-Ethene """
TAGL['%s-%s-monoB-unCP' % (dbse, '30'                    )] = """Monomer B from Benzene-Ethene """
TAGL['%s-%s'            % (dbse, '31'                    )] = """Uracil-Ethene """
TAGL['%s-%s-dimer'      % (dbse, '31'                    )] = """Dimer from Uracil-Ethene """
TAGL['%s-%s-monoA-CP'   % (dbse, '31'                    )] = """Monomer A from Uracil-Ethene """
TAGL['%s-%s-monoB-CP'   % (dbse, '31'                    )] = """Monomer B from Uracil-Ethene """
TAGL['%s-%s-monoA-unCP' % (dbse, '31'                    )] = """Monomer A from Uracil-Ethene """
TAGL['%s-%s-monoB-unCP' % (dbse, '31'                    )] = """Monomer B from Uracil-Ethene """
TAGL['%s-%s'            % (dbse, '32'                    )] = """Uracil-Ethyne """
TAGL['%s-%s-dimer'      % (dbse, '32'                    )] = """Dimer from Uracil-Ethyne """
TAGL['%s-%s-monoA-CP'   % (dbse, '32'                    )] = """Monomer A from Uracil-Ethyne """
TAGL['%s-%s-monoB-CP'   % (dbse, '32'                    )] = """Monomer B from Uracil-Ethyne """
TAGL['%s-%s-monoA-unCP' % (dbse, '32'                    )] = """Monomer A from Uracil-Ethyne """
TAGL['%s-%s-monoB-unCP' % (dbse, '32'                    )] = """Monomer B from Uracil-Ethyne """
TAGL['%s-%s'            % (dbse, '33'                    )] = """Pyridine-Ethene """
TAGL['%s-%s-dimer'      % (dbse, '33'                    )] = """Dimer from Pyridine-Ethene """
TAGL['%s-%s-monoA-CP'   % (dbse, '33'                    )] = """Monomer A from Pyridine-Ethene """
TAGL['%s-%s-monoB-CP'   % (dbse, '33'                    )] = """Monomer B from Pyridine-Ethene """
TAGL['%s-%s-monoA-unCP' % (dbse, '33'                    )] = """Monomer A from Pyridine-Ethene """
TAGL['%s-%s-monoB-unCP' % (dbse, '33'                    )] = """Monomer B from Pyridine-Ethene """
TAGL['%s-%s'            % (dbse, '34'                    )] = """Pentane Dimer """
TAGL['%s-%s-dimer'      % (dbse, '34'                    )] = """Dimer from Pentane Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '34'                    )] = """Monomer A from Pentane Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '34'                    )] = """Monomer B from Pentane Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '34'                    )] = """Monomer A from Pentane Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '34'                    )] = """Monomer B from Pentane Dimer """
TAGL['%s-%s'            % (dbse, '35'                    )] = """Neopentane-Pentane """
TAGL['%s-%s-dimer'      % (dbse, '35'                    )] = """Dimer from Neopentane-Pentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '35'                    )] = """Monomer A from Neopentane-Pentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '35'                    )] = """Monomer B from Neopentane-Pentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '35'                    )] = """Monomer A from Neopentane-Pentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '35'                    )] = """Monomer B from Neopentane-Pentane """
TAGL['%s-%s'            % (dbse, '36'                    )] = """Neopentane Dimer """
TAGL['%s-%s-dimer'      % (dbse, '36'                    )] = """Dimer from Neopentane Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '36'                    )] = """Monomer A from Neopentane Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '36'                    )] = """Monomer B from Neopentane Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '36'                    )] = """Monomer A from Neopentane Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '36'                    )] = """Monomer B from Neopentane Dimer """
TAGL['%s-%s'            % (dbse, '37'                    )] = """Cyclopentane-Neopentane """
TAGL['%s-%s-dimer'      % (dbse, '37'                    )] = """Dimer from Cyclopentane-Neopentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '37'                    )] = """Monomer A from Cyclopentane-Neopentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '37'                    )] = """Monomer B from Cyclopentane-Neopentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '37'                    )] = """Monomer A from Cyclopentane-Neopentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '37'                    )] = """Monomer B from Cyclopentane-Neopentane """
TAGL['%s-%s'            % (dbse, '38'                    )] = """Cyclopentane Dimer """
TAGL['%s-%s-dimer'      % (dbse, '38'                    )] = """Dimer from Cyclopentane Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '38'                    )] = """Monomer A from Cyclopentane Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '38'                    )] = """Monomer B from Cyclopentane Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '38'                    )] = """Monomer A from Cyclopentane Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '38'                    )] = """Monomer B from Cyclopentane Dimer """
TAGL['%s-%s'            % (dbse, '39'                    )] = """Benzene-Cyclopentane """
TAGL['%s-%s-dimer'      % (dbse, '39'                    )] = """Dimer from Benzene-Cyclopentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '39'                    )] = """Monomer A from Benzene-Cyclopentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '39'                    )] = """Monomer B from Benzene-Cyclopentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '39'                    )] = """Monomer A from Benzene-Cyclopentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '39'                    )] = """Monomer B from Benzene-Cyclopentane """
TAGL['%s-%s'            % (dbse, '40'                    )] = """Benzene-Neopentane """
TAGL['%s-%s-dimer'      % (dbse, '40'                    )] = """Dimer from Benzene-Neopentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '40'                    )] = """Monomer A from Benzene-Neopentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '40'                    )] = """Monomer B from Benzene-Neopentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '40'                    )] = """Monomer A from Benzene-Neopentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '40'                    )] = """Monomer B from Benzene-Neopentane """
TAGL['%s-%s'            % (dbse, '41'                    )] = """Uracil-Pentane """
TAGL['%s-%s-dimer'      % (dbse, '41'                    )] = """Dimer from Uracil-Pentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '41'                    )] = """Monomer A from Uracil-Pentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '41'                    )] = """Monomer B from Uracil-Pentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '41'                    )] = """Monomer A from Uracil-Pentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '41'                    )] = """Monomer B from Uracil-Pentane """
TAGL['%s-%s'            % (dbse, '42'                    )] = """Uracil-Cyclopentane """
TAGL['%s-%s-dimer'      % (dbse, '42'                    )] = """Dimer from Uracil-Cyclopentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '42'                    )] = """Monomer A from Uracil-Cyclopentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '42'                    )] = """Monomer B from Uracil-Cyclopentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '42'                    )] = """Monomer A from Uracil-Cyclopentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '42'                    )] = """Monomer B from Uracil-Cyclopentane """
TAGL['%s-%s'            % (dbse, '43'                    )] = """Uracil-Neopentane """
TAGL['%s-%s-dimer'      % (dbse, '43'                    )] = """Dimer from Uracil-Neopentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '43'                    )] = """Monomer A from Uracil-Neopentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '43'                    )] = """Monomer B from Uracil-Neopentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '43'                    )] = """Monomer A from Uracil-Neopentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '43'                    )] = """Monomer B from Uracil-Neopentane """
TAGL['%s-%s'            % (dbse, '44'                    )] = """Ethene-Pentane """
TAGL['%s-%s-dimer'      % (dbse, '44'                    )] = """Dimer from Ethene-Pentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '44'                    )] = """Monomer A from Ethene-Pentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '44'                    )] = """Monomer B from Ethene-Pentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '44'                    )] = """Monomer A from Ethene-Pentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '44'                    )] = """Monomer B from Ethene-Pentane """
TAGL['%s-%s'            % (dbse, '45'                    )] = """Ethyne-Pentane """
TAGL['%s-%s-dimer'      % (dbse, '45'                    )] = """Dimer from Ethyne-Pentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '45'                    )] = """Monomer A from Ethyne-Pentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '45'                    )] = """Monomer B from Ethyne-Pentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '45'                    )] = """Monomer A from Ethyne-Pentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '45'                    )] = """Monomer B from Ethyne-Pentane """
TAGL['%s-%s'            % (dbse, '46'                    )] = """N-methylacetamide-Pentane """
TAGL['%s-%s-dimer'      % (dbse, '46'                    )] = """Dimer from N-methylacetamide-Pentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '46'                    )] = """Monomer A from N-methylacetamide-Pentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '46'                    )] = """Monomer B from N-methylacetamide-Pentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '46'                    )] = """Monomer A from N-methylacetamide-Pentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '46'                    )] = """Monomer B from N-methylacetamide-Pentane """
TAGL['%s-%s'            % (dbse, '47'                    )] = """Benzene Dimer, CH-pi """
TAGL['%s-%s-dimer'      % (dbse, '47'                    )] = """Dimer from Benzene Dimer, CH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '47'                    )] = """Monomer A from Benzene Dimer, CH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '47'                    )] = """Monomer B from Benzene Dimer, CH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '47'                    )] = """Monomer A from Benzene Dimer, CH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '47'                    )] = """Monomer B from Benzene Dimer, CH-pi """
TAGL['%s-%s'            % (dbse, '48'                    )] = """Pyridine Dimer, CH-pi """
TAGL['%s-%s-dimer'      % (dbse, '48'                    )] = """Dimer from Pyridine Dimer, CH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '48'                    )] = """Monomer A from Pyridine Dimer, CH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '48'                    )] = """Monomer B from Pyridine Dimer, CH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '48'                    )] = """Monomer A from Pyridine Dimer, CH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '48'                    )] = """Monomer B from Pyridine Dimer, CH-pi """
TAGL['%s-%s'            % (dbse, '49'                    )] = """Benzene-Pyridine, CH-pi """
TAGL['%s-%s-dimer'      % (dbse, '49'                    )] = """Dimer from Benzene-Pyridine, CH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '49'                    )] = """Monomer A from Benzene-Pyridine, CH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '49'                    )] = """Monomer B from Benzene-Pyridine, CH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '49'                    )] = """Monomer A from Benzene-Pyridine, CH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '49'                    )] = """Monomer B from Benzene-Pyridine, CH-pi """
TAGL['%s-%s'            % (dbse, '50'                    )] = """Benzene-Ethyne, CH-pi """
TAGL['%s-%s-dimer'      % (dbse, '50'                    )] = """Dimer from Benzene-Ethyne, CH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '50'                    )] = """Monomer A from Benzene-Ethyne, CH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '50'                    )] = """Monomer B from Benzene-Ethyne, CH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '50'                    )] = """Monomer A from Benzene-Ethyne, CH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '50'                    )] = """Monomer B from Benzene-Ethyne, CH-pi """
TAGL['%s-%s'            % (dbse, '51'                    )] = """Ethyne Dimer, CH-pi """
TAGL['%s-%s-dimer'      % (dbse, '51'                    )] = """Dimer from Ethyne Dimer, CH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '51'                    )] = """Monomer A from Ethyne Dimer, CH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '51'                    )] = """Monomer B from Ethyne Dimer, CH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '51'                    )] = """Monomer A from Ethyne Dimer, CH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '51'                    )] = """Monomer B from Ethyne Dimer, CH-pi """
TAGL['%s-%s'            % (dbse, '52'                    )] = """Benzene-Acetic Acid, OH-pi """
TAGL['%s-%s-dimer'      % (dbse, '52'                    )] = """Dimer from Benzene-Acetic Acid, OH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '52'                    )] = """Monomer A from Benzene-Acetic Acid, OH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '52'                    )] = """Monomer B from Benzene-Acetic Acid, OH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '52'                    )] = """Monomer A from Benzene-Acetic Acid, OH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '52'                    )] = """Monomer B from Benzene-Acetic Acid, OH-pi """
TAGL['%s-%s'            % (dbse, '53'                    )] = """Benzene-Acetamide, NH-pi """
TAGL['%s-%s-dimer'      % (dbse, '53'                    )] = """Dimer from Benzene-Acetamide, NH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '53'                    )] = """Monomer A from Benzene-Acetamide, NH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '53'                    )] = """Monomer B from Benzene-Acetamide, NH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '53'                    )] = """Monomer A from Benzene-Acetamide, NH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '53'                    )] = """Monomer B from Benzene-Acetamide, NH-pi """
TAGL['%s-%s'            % (dbse, '54'                    )] = """Benzene-Water, OH-pi """
TAGL['%s-%s-dimer'      % (dbse, '54'                    )] = """Dimer from Benzene-Water, OH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '54'                    )] = """Monomer A from Benzene-Water, OH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '54'                    )] = """Monomer B from Benzene-Water, OH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '54'                    )] = """Monomer A from Benzene-Water, OH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '54'                    )] = """Monomer B from Benzene-Water, OH-pi """
TAGL['%s-%s'            % (dbse, '55'                    )] = """Benzene-Methanol, OH-pi """
TAGL['%s-%s-dimer'      % (dbse, '55'                    )] = """Dimer from Benzene-Methanol, OH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '55'                    )] = """Monomer A from Benzene-Methanol, OH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '55'                    )] = """Monomer B from Benzene-Methanol, OH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '55'                    )] = """Monomer A from Benzene-Methanol, OH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '55'                    )] = """Monomer B from Benzene-Methanol, OH-pi """
TAGL['%s-%s'            % (dbse, '56'                    )] = """Benzene-Methylamine, NH-pi """
TAGL['%s-%s-dimer'      % (dbse, '56'                    )] = """Dimer from Benzene-Methylamine, NH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '56'                    )] = """Monomer A from Benzene-Methylamine, NH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '56'                    )] = """Monomer B from Benzene-Methylamine, NH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '56'                    )] = """Monomer A from Benzene-Methylamine, NH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '56'                    )] = """Monomer B from Benzene-Methylamine, NH-pi """
TAGL['%s-%s'            % (dbse, '57'                    )] = """Benzene-N-methylacetamide, NH-pi """
TAGL['%s-%s-dimer'      % (dbse, '57'                    )] = """Dimer from Benzene-N-methylacetamide, NH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '57'                    )] = """Monomer A from Benzene-N-methylacetamide, NH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '57'                    )] = """Monomer B from Benzene-N-methylacetamide, NH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '57'                    )] = """Monomer A from Benzene-N-methylacetamide, NH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '57'                    )] = """Monomer B from Benzene-N-methylacetamide, NH-pi """
TAGL['%s-%s'            % (dbse, '58'                    )] = """Pyridine Dimer, CH-N """
TAGL['%s-%s-dimer'      % (dbse, '58'                    )] = """Dimer from Pyridine Dimer, CH-N """
TAGL['%s-%s-monoA-CP'   % (dbse, '58'                    )] = """Monomer A from Pyridine Dimer, CH-N """
TAGL['%s-%s-monoB-CP'   % (dbse, '58'                    )] = """Monomer B from Pyridine Dimer, CH-N """
TAGL['%s-%s-monoA-unCP' % (dbse, '58'                    )] = """Monomer A from Pyridine Dimer, CH-N """
TAGL['%s-%s-monoB-unCP' % (dbse, '58'                    )] = """Monomer B from Pyridine Dimer, CH-N """
TAGL['%s-%s'            % (dbse, '59'                    )] = """Ethyne-Water, CH-O """
TAGL['%s-%s-dimer'      % (dbse, '59'                    )] = """Dimer from Ethyne-Water, CH-O """
TAGL['%s-%s-monoA-CP'   % (dbse, '59'                    )] = """Monomer A from Ethyne-Water, CH-O """
TAGL['%s-%s-monoB-CP'   % (dbse, '59'                    )] = """Monomer B from Ethyne-Water, CH-O """
TAGL['%s-%s-monoA-unCP' % (dbse, '59'                    )] = """Monomer A from Ethyne-Water, CH-O """
TAGL['%s-%s-monoB-unCP' % (dbse, '59'                    )] = """Monomer B from Ethyne-Water, CH-O """
TAGL['%s-%s'            % (dbse, '60'                    )] = """Ethyne-Acetic Acid, OH-pi """
TAGL['%s-%s-dimer'      % (dbse, '60'                    )] = """Dimer from Ethyne-Acetic Acid, OH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '60'                    )] = """Monomer A from Ethyne-Acetic Acid, OH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '60'                    )] = """Monomer B from Ethyne-Acetic Acid, OH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '60'                    )] = """Monomer A from Ethyne-Acetic Acid, OH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '60'                    )] = """Monomer B from Ethyne-Acetic Acid, OH-pi """
TAGL['%s-%s'            % (dbse, '61'                    )] = """Pentane-Acetic Acid """
TAGL['%s-%s-dimer'      % (dbse, '61'                    )] = """Dimer from Pentane-Acetic Acid """
TAGL['%s-%s-monoA-CP'   % (dbse, '61'                    )] = """Monomer A from Pentane-Acetic Acid """
TAGL['%s-%s-monoB-CP'   % (dbse, '61'                    )] = """Monomer B from Pentane-Acetic Acid """
TAGL['%s-%s-monoA-unCP' % (dbse, '61'                    )] = """Monomer A from Pentane-Acetic Acid """
TAGL['%s-%s-monoB-unCP' % (dbse, '61'                    )] = """Monomer B from Pentane-Acetic Acid """
TAGL['%s-%s'            % (dbse, '62'                    )] = """Pentane-Acetamide """
TAGL['%s-%s-dimer'      % (dbse, '62'                    )] = """Dimer from Pentane-Acetamide """
TAGL['%s-%s-monoA-CP'   % (dbse, '62'                    )] = """Monomer A from Pentane-Acetamide """
TAGL['%s-%s-monoB-CP'   % (dbse, '62'                    )] = """Monomer B from Pentane-Acetamide """
TAGL['%s-%s-monoA-unCP' % (dbse, '62'                    )] = """Monomer A from Pentane-Acetamide """
TAGL['%s-%s-monoB-unCP' % (dbse, '62'                    )] = """Monomer B from Pentane-Acetamide """
TAGL['%s-%s'            % (dbse, '63'                    )] = """Benzene-Acetic Acid """
TAGL['%s-%s-dimer'      % (dbse, '63'                    )] = """Dimer from Benzene-Acetic Acid """
TAGL['%s-%s-monoA-CP'   % (dbse, '63'                    )] = """Monomer A from Benzene-Acetic Acid """
TAGL['%s-%s-monoB-CP'   % (dbse, '63'                    )] = """Monomer B from Benzene-Acetic Acid """
TAGL['%s-%s-monoA-unCP' % (dbse, '63'                    )] = """Monomer A from Benzene-Acetic Acid """
TAGL['%s-%s-monoB-unCP' % (dbse, '63'                    )] = """Monomer B from Benzene-Acetic Acid """
TAGL['%s-%s'            % (dbse, '64'                    )] = """N-methylacetamide-Ethene """
TAGL['%s-%s-dimer'      % (dbse, '64'                    )] = """Dimer from N-methylacetamide-Ethene """
TAGL['%s-%s-monoA-CP'   % (dbse, '64'                    )] = """Monomer A from N-methylacetamide-Ethene """
TAGL['%s-%s-monoB-CP'   % (dbse, '64'                    )] = """Monomer B from N-methylacetamide-Ethene """
TAGL['%s-%s-monoA-unCP' % (dbse, '64'                    )] = """Monomer A from N-methylacetamide-Ethene """
TAGL['%s-%s-monoB-unCP' % (dbse, '64'                    )] = """Monomer B from N-methylacetamide-Ethene """
TAGL['%s-%s'            % (dbse, '65'                    )] = """Pyridine-Ethyne """
TAGL['%s-%s-dimer'      % (dbse, '65'                    )] = """Dimer from Pyridine-Ethyne """
TAGL['%s-%s-monoA-CP'   % (dbse, '65'                    )] = """Monomer A from Pyridine-Ethyne """
TAGL['%s-%s-monoB-CP'   % (dbse, '65'                    )] = """Monomer B from Pyridine-Ethyne """
TAGL['%s-%s-monoA-unCP' % (dbse, '65'                    )] = """Monomer A from Pyridine-Ethyne """
TAGL['%s-%s-monoB-unCP' % (dbse, '65'                    )] = """Monomer B from Pyridine-Ethyne """
TAGL['%s-%s'            % (dbse, '66'                    )] = """Methylamine-Pyridine """
TAGL['%s-%s-dimer'      % (dbse, '66'                    )] = """Dimer from Methylamine-Pyridine """
TAGL['%s-%s-monoA-CP'   % (dbse, '66'                    )] = """Monomer A from Methylamine-Pyridine """
TAGL['%s-%s-monoB-CP'   % (dbse, '66'                    )] = """Monomer B from Methylamine-Pyridine """
TAGL['%s-%s-monoA-unCP' % (dbse, '66'                    )] = """Monomer A from Methylamine-Pyridine """
TAGL['%s-%s-monoB-unCP' % (dbse, '66'                    )] = """Monomer B from Methylamine-Pyridine """

# <<< Molecule Specifications >>>
monoA_unCP = 'monoA = dimer.extract_subsets(1)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_unCP = 'monoB = dimer.extract_subsets(2)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'
monoA_CP   = 'monoA = dimer.extract_subsets(1,2)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_CP   = 'monoB = dimer.extract_subsets(2,1)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'

S66_1 = input.process_input("""
molecule dimer {
0 1
O       -0.70219605    -0.05606026     0.00994226
H       -1.02219322     0.84677578    -0.01148871
H        0.25752106     0.04212150     0.00521900
--
0 1
O        2.22087107     0.02671679     0.00062048
H        2.59749268    -0.41166327     0.76674486
H        2.59313538    -0.44949618    -0.74478203
units angstrom
}
""", 0)

S66_2 = input.process_input("""
molecule dimer {
0 1
O       -0.52532979    -0.05097108    -0.31451686
H       -0.94200663     0.74790163     0.01125282
H        0.40369652     0.05978598    -0.07356837
--
0 1
O        2.31663329     0.04550085     0.07185839
H        2.68461611    -0.52657655     0.74938672
C        2.78163836    -0.42612907    -1.19030072
H        2.35082127     0.22496462    -1.94341475
H        3.86760205    -0.37533621    -1.26461265
H        2.45329574    -1.44599856    -1.38938136
units angstrom
}
""", 0)

S66_3 = input.process_input("""
molecule dimer {
0 1
O       -0.68746490    -0.11174433    -0.01962547
H       -1.04612154     0.77593821     0.01270684
H        0.27404252     0.02585065    -0.00349726
--
0 1
N        2.23397617     0.10318260     0.00585368
H        2.52934060    -0.44945538    -0.78893718
H        2.54405666    -0.40753849     0.82271317
C        2.89331145     1.41154656    -0.03438796
H        2.58276902     1.99327152     0.83012746
H        3.98462074     1.37225159    -0.04334363
H        2.56659917     1.94746403    -0.92221177
units angstrom
}
""", 0)

S66_4 = input.process_input("""
molecule dimer {
0 1
O       -0.39201845    -0.38471874     0.07607132
H       -0.91146085     0.41381204     0.17764877
H        0.52490382    -0.06848469     0.09051136
--
0 1
C        2.19770521    -2.24540349    -0.23031325
H        2.84766805    -3.10651537    -0.36322864
H        1.51672924    -2.16793143    -1.07417853
H        1.58468831    -2.38419948     0.65669511
C        2.95243729    -0.94739061    -0.09771974
O        2.37572184     0.12790424     0.05886900
N        4.30307041    -1.04489330    -0.16233771
H        4.70402204    -1.95542728    -0.29185281
C        5.17131253     0.10707716    -0.05289463
H        4.53481840     0.97537761     0.08188998
H        5.83690203     0.01562196     0.80319825
H        5.76577825     0.23649765    -0.95515382
units angstrom
}
""", 0)

S66_5 = input.process_input("""
molecule dimer {
0 1
O       -0.63613493    -0.02328241     0.28059932
H        0.30809737    -0.04707875     0.07646369
C       -1.15206541    -1.31128778     0.01525955
H       -2.20994502    -1.29626539     0.26395586
H       -1.05661024    -1.59267086    -1.03619061
H       -0.67483575    -2.08627276     0.62051145
--
0 1
O        2.21041928    -0.12212177    -0.01210270
H        2.67920859     0.49226275    -0.58176865
C        2.71925320     0.03489717     1.30961462
H        2.16568412    -0.65329926     1.93974550
H        3.77824931    -0.21554173     1.36633776
H        2.56681356     1.04559122     1.68750717
units angstrom
}
""", 0)

S66_6 = input.process_input("""
molecule dimer {
0 1
O       -0.70692019     0.04583037     0.00638610
H        0.26562361     0.07171014     0.00133929
C       -1.07667067    -1.31391581     0.00161428
H       -2.16292358    -1.36319577     0.00586542
H       -0.72340594    -1.84465168    -0.88774350
H       -0.71607978    -1.85282083     0.88307978
--
0 1
N        2.20127244    -0.03642087    -0.00333839
H        2.57189199     0.47135563     0.78979400
H        2.57201528     0.42791769    -0.82259722
C        2.67902438    -1.42245432     0.03412282
H        2.28713954    -1.95647960    -0.82806891
H        3.76573553    -1.52918949     0.03715731
H        2.28689798    -1.90918449     0.92375496
units angstrom
}
""", 0)

S66_7 = input.process_input("""
molecule dimer {
0 1
O       -0.20877739    -0.21687067    -1.03240597
H        0.71112593    -0.38689175    -0.77396240
C       -1.02217337    -0.74117114    -0.00545419
H       -2.05749119    -0.53870733    -0.26859725
H       -0.90774336    -1.82182632     0.10853710
H       -0.82463111    -0.27549472     0.96464547
--
0 1
C        1.97349049     1.90322403     0.43230118
H        2.47988412     2.86467311     0.39743082
H        1.56294637     1.75708815     1.43017782
H        1.14384269     1.89371075    -0.26920435
C        2.88912087     0.74828521     0.11638497
O        2.46492608    -0.37162558    -0.16869657
N        4.21525779     1.01000949     0.17558433
H        4.51327024     1.92043762     0.47327152
C        5.19766382    -0.03010182    -0.04715949
H        4.84110663    -0.68103914    -0.83933645
H        6.13803306     0.42342202    -0.34567319
H        5.35717393    -0.63462872     0.84491605
units angstrom
}
""", 0)

S66_8 = input.process_input("""
molecule dimer {
0 1
O       -0.78656202     0.04516844    -0.00718912
H        0.17770677     0.01269590    -0.00683539
C       -1.24799094    -1.29028354     0.00108362
H       -2.33427744    -1.25889710     0.00022120
H       -0.92596575    -1.84976810    -0.88044538
H       -0.92702783    -1.83846288     0.89007652
--
0 1
O        2.12888314    -0.05133660    -0.00474093
H        2.56808728     0.33681560    -0.76461362
H        2.56676744     0.35126768     0.74834860
units angstrom
}
""", 0)

S66_9 = input.process_input("""
molecule dimer {
0 1
N       -0.89345122    -0.04384432    -0.04299745
H        0.09694826    -0.25605945    -0.07106993
H       -1.36843879    -0.93339065     0.03383773
C       -1.17578248     0.75790769     1.14523719
H       -2.24162660     0.97221601     1.19502464
H       -0.88078955     0.30424674     2.09720910
H       -0.66300572     1.71432940     1.06080916
--
0 1
O        2.28445953    -0.04747650     0.02782522
H        2.56648565     0.32247227    -0.81203886
C        2.67037338     0.86410776     1.04726138
H        2.34719033     0.43447509     1.99032792
H        3.75142862     1.00319123     1.08630135
H        2.19189882     1.83770561     0.93208484
units angstrom
}
""", 0)

S66_10 = input.process_input("""
molecule dimer {
0 1
N       -0.63864138     0.47091637     0.04456848
H        0.18995436    -0.11393716    -0.00577361
H       -1.30046894     0.08125680    -0.61366848
C       -1.19865882     0.39139858     1.39194660
H       -2.09273777     1.00924471     1.45316749
H       -1.46274551    -0.61584367     1.72945219
H       -0.48027554     0.79867491     2.10108731
--
0 1
N        2.39889347    -0.45552115     0.19704452
H        2.69516214    -0.18098342    -0.73094072
H        3.02244314    -1.20321147     0.47223938
C        2.55912345     0.67968944     1.11071982
H        2.28893315     0.36499366     2.11637293
H        3.56653376     1.10146600     1.14769156
H        1.86658307     1.46546492     0.81806258
units angstrom
}
""", 0)

S66_11 = input.process_input("""
molecule dimer {
0 1
N       -0.56970824     0.81437245     0.10109775
H        0.13087774     0.56141065    -0.58761455
H       -1.46125215     0.52691480    -0.28042996
C       -0.30551437     0.06571030     1.32879173
H       -1.05714948     0.31427017     2.07595940
H       -0.28802353    -1.02229248     1.21484626
H        0.66045772     0.36850913     1.73024224
--
0 1
C        2.25689155     2.69009990    -0.14932730
H        2.38151002     3.10127663    -1.14837163
H        2.76346292     3.33109245     0.56845722
H        1.19047979     2.66357037     0.06909413
C        2.76888324     1.27230222    -0.14703327
O        2.30890335     0.40656580    -0.88620788
N        3.75536621     0.99926987     0.74529744
H        4.15512723     1.75420265     1.27065019
C        4.34381155    -0.32032067     0.82279701
H        3.55563493    -1.06165082     0.72977641
H        5.06507133    -0.49231605     0.02425262
H        4.83846506    -0.43618886     1.78273654
units angstrom
}
""", 0)

S66_12 = input.process_input("""
molecule dimer {
0 1
N       -0.53346397    -0.27959351     0.10699576
H       -0.62915138    -1.24842455     0.38284867
H       -1.12260363    -0.16615944    -0.70776410
C       -1.01690943     0.58848610     1.18737346
H       -0.91275967     1.62555174     0.87952116
H       -2.05473726     0.41508213     1.47850360
H       -0.38502338     0.44880090     2.06061419
--
0 1
O        2.09326841     0.91731136     0.21209725
H        1.27575101     0.42103887     0.03894435
H        2.67516986     0.65881349    -0.50364884
units angstrom
}
""", 0)

S66_13 = input.process_input("""
molecule dimer {
0 1
C       -0.84931672    -0.33949876     2.49171664
H        0.18434396    -0.01104732     2.41618542
H       -0.88249791    -1.34205140     2.91270310
H       -1.39080263     0.31687828     3.16842897
C       -1.56403192    -0.35332311     1.15947545
O       -2.74952638    -0.65153776     1.05676087
N       -0.80165352    -0.02735461     0.08834167
H        0.16118756     0.24036035     0.21871364
C       -1.38534986    -0.00235149    -1.23413683
H       -1.89161720    -0.94280123    -1.44009631
H       -2.11997230     0.79621180    -1.33087952
H       -0.59464593     0.14957065    -1.96312772
--
0 1
O        2.13706570     0.25201737     0.45371880
H        2.85792051     0.87931700     0.54413361
C        2.65614986    -1.05334828     0.68760059
H        1.82357836    -1.74213597     0.58202402
H        3.42228862    -1.32234103    -0.03928018
H        3.06424691    -1.15479748     1.69323508
units angstrom
}
""", 0)

S66_14 = input.process_input("""
molecule dimer {
0 1
C       -0.77857334    -0.46332064     2.49038768
H        0.22474462    -0.05095294     2.41348355
H       -0.72247994    -1.48709180     2.85458464
H       -1.35190757     0.11081693     3.21368365
C       -1.52050259    -0.45662769     1.17232500
O       -2.70083521    -0.78358573     1.08959682
N       -0.79195361    -0.06964048     0.10058937
H        0.19411165     0.14570790     0.20292464
C       -1.39779834    -0.05608245    -1.21131793
H       -2.31492801     0.52889121    -1.19970991
H       -0.69880422     0.38726130    -1.91536621
H       -1.65298232    -1.06152895    -1.54543495
--
0 1
N        2.23828822     0.25457428     0.28251924
H        2.64195454     0.79449381     1.03771933
H        2.65629209     0.62195553    -0.56312668
C        2.61059106    -1.15660854     0.43627199
H        2.18430366    -1.72764112    -0.38510346
H        3.68598970    -1.34329798     0.46205539
H        2.17611849    -1.54101555     1.35610799
units angstrom
}
""", 0)

S66_15 = input.process_input("""
molecule dimer {
0 1
C       -0.70150294    -0.29062770     2.40688440
H       -1.18329596     0.39564777     3.09887422
H        0.34956157    -0.03032157     2.30783303
H       -0.79405685    -1.29160545     2.82403929
C       -1.44854625    -0.24487664     1.09181530
O       -2.66045000    -0.42847909     1.03434577
N       -0.67005656     0.00591656     0.00977691
H        0.32667532     0.12256396     0.14159284
C       -1.22705457     0.08979374    -1.31996754
H       -2.29202426    -0.10650119    -1.24087756
H       -1.07780169     1.07994030    -1.74854354
H       -0.77662849    -0.64799919    -1.98337273
--
0 1
C        2.04177491    -2.35169797     0.68639761
H        2.59999972    -3.26170120     0.48048961
H        1.11308306    -2.35822742     0.12207220
H        1.78255599    -2.32825127     1.74333861
C        2.80941086    -1.09728593     0.35016088
O        2.26422421     0.00415088     0.29318848
N        4.13616907    -1.26609970     0.13641291
H        4.51249037    -2.19334539     0.21317023
C        5.02340725    -0.15963372    -0.15253563
H        4.40921487     0.73117605    -0.23235934
H        5.75082180    -0.02016799     0.64486768
H        5.54839755    -0.31961545    -1.09167796
units angstrom
}
""", 0)

S66_16 = input.process_input("""
molecule dimer {
0 1
C       -0.72430464    -0.70493582     2.28386786
H        0.33531828    -0.62994325     2.05318235
H       -0.95169666    -1.71198961     2.62565146
H       -0.96962784    -0.02207955     3.09376537
C       -1.61493501    -0.38742925     1.10406897
O       -2.83732387    -0.41502209     1.19413277
N       -0.95342037    -0.07640442    -0.04081980
H        0.05380860    -0.07556651    -0.03664022
C       -1.65812397     0.25009358    -1.25855306
H       -2.72037197     0.17694444    -1.04665270
H       -1.43030493     1.26296263    -1.58809384
H       -1.40562611    -0.44433518    -2.05858358
--
0 1
O        2.10277707    -0.05840697    -0.15507669
H        2.66775436    -0.77136560    -0.46027609
H        2.68252869     0.70578659    -0.13117819
units angstrom
}
""", 0)

S66_17 = input.process_input("""
molecule dimer {
0 1
N       -0.72999913     0.02276763     0.00091465
H        0.29842255     0.07400447     0.00162304
C       -1.29682453    -1.24042682     0.00150234
O       -0.59409886    -2.25351751     0.00263371
C       -2.74362229    -1.26233170     0.00047938
H       -3.24959045    -2.21183517     0.00083311
C       -3.42201997    -0.09590921    -0.00092259
H       -4.50089709    -0.04921603    -0.00174546
N       -2.77483684     1.10540895    -0.00141807
H       -3.28383807     1.97387739    -0.00248574
C       -1.39147866     1.23701978    -0.00052538
O       -0.83984371     2.31703528    -0.00100125
--
0 1
N        4.14382946    -1.08570382     0.00049928
H        4.59107325    -0.17913062     0.00088609
C        4.99987723    -2.20032161    -0.00100060
O        6.20932926    -2.04861719    -0.00174980
C        4.28565880    -3.46249515    -0.00150500
H        4.85224335    -4.37752590    -0.00264363
C        2.93548983    -3.46631302    -0.00054490
H        2.35852659    -4.37927779    -0.00086358
N        2.19749842    -2.31543218     0.00090551
H        1.17116216    -2.33687498     0.00158258
C        2.77026935    -1.07076714     0.00145616
O        2.11994847    -0.02954883     0.00269255
units angstrom
}
""", 0)

S66_18 = input.process_input("""
molecule dimer {
0 1
O       -0.55283102    -0.10169749    -0.00049879
H       -0.87175963     0.80179220     0.00014440
H        0.41265950    -0.00183225    -0.00025181
--
0 1
N        2.36402099     0.09662268     0.00014680
C        3.05992763     0.06265189     1.14489465
H        2.47525508     0.08626283     2.05576267
C        4.44895122    -0.00253054     1.19489071
H        4.95485760    -0.02738470     2.14921983
C        5.16011436    -0.03565634    -0.00002044
H        6.23995431    -0.08742989    -0.00010086
C        4.44880607    -0.00259720    -1.19482173
H        4.95460301    -0.02747022    -2.14922033
C        3.05977605     0.06259779    -1.14467547
H        2.47500717     0.08619845    -2.05546803
units angstrom
}
""", 0)

S66_19 = input.process_input("""
molecule dimer {
0 1
O       -0.62765177     0.08746727     0.00147128
H        0.34360203     0.12230333    -0.00060045
C       -0.97793123    -1.27855601     0.00123841
H       -2.06339209    -1.34204332     0.00500898
H       -0.61488369    -1.80637584    -0.88538395
H       -0.60864033    -1.80823682     0.88417273
--
0 1
N        2.27233665     0.01643230    -0.00162684
C        2.96870504    -0.00800303    -1.14634644
H        2.38422645     0.01522051    -2.05732188
C        4.35834211    -0.05774589    -1.19503169
H        4.86569445    -0.07503793    -2.14881442
C        5.06871533    -0.08345851     0.00058133
H        6.14905134    -0.12122326     0.00143063
C        4.35646788    -0.05843740     1.19512119
H        4.86226662    -0.07626173     2.14960688
C        2.96691424    -0.00868772     1.14416710
H        2.38090845     0.01398671     2.05428579
units angstrom
}
""", 0)

S66_20 = input.process_input("""
molecule dimer {
0 1
C       -1.06170920     1.29714057     0.29206000
O       -0.35816112     2.27045861     0.53181267
O       -0.58930352     0.09491776     0.00378881
H        0.40443566     0.12772262     0.01841184
C       -2.55842780     1.34254982     0.29625732
H       -2.89599798     2.34746400     0.51831634
H       -2.93288928     1.02239045    -0.67299555
H       -2.93721196     0.64491043     1.03955708
--
0 1
C        2.78934845     1.10841924     0.27118376
O        2.08573008     0.13510475     0.03139616
O        2.31692211     2.31085463     0.55896223
H        1.32313357     2.27795640     0.54456172
C        4.28606090     1.06251650     0.26921936
H        4.62364046     0.06119730     0.03169387
H        4.66755944     1.77286944    -0.46024953
H        4.65757721     1.36521101     1.24527472
units angstrom
}
""", 0)

S66_21 = input.process_input("""
molecule dimer {
0 1
C       -1.30974974     1.18017617    -0.02517034
O       -0.72530044     2.15514767     0.45271335
N       -0.66562116     0.09505470    -0.49199449
H        0.35458266     0.05144817    -0.45930922
H       -1.18362704    -0.67359969    -0.87075610
C       -2.81671934     1.15599865    -0.11060597
H       -3.22062895     1.26254146     0.89308239
H       -3.20942754     0.24863402    -0.56190009
H       -3.14315813     2.01659563    -0.68889311
--
0 1
C        2.77960183     1.06388568     0.13435724
O        2.19518007     0.08986525    -0.34537373
N        2.13551426     2.14862891     0.60220379
H        1.11540890     2.19306669     0.56790248
H        2.65353833     2.91659011     0.98232444
C        4.28660101     1.08817006     0.21958232
H        4.67847207     1.98781958     0.68676633
H        4.69015720     1.00062503    -0.78619798
H        4.61437977     0.21759516     0.78176266
units angstrom
}
""", 0)

S66_22 = input.process_input("""
molecule dimer {
0 1
C       -1.11362611     1.32702009     0.27516705
O       -0.46708264     2.34938778     0.46153746
O       -0.57808939     0.13692049     0.04961747
H        0.41332036     0.20325661     0.05548711
C       -2.61142469     1.28618957     0.27736131
H       -3.00664872     2.27688545     0.46578983
H       -2.96425623     0.91525868    -0.68200123
H       -2.95311421     0.59179821     1.04124041
--
0 1
N        4.18869738     1.08795338     0.18288157
H        4.58190249     0.17256315     0.01116215
C        5.11022529     2.13606900     0.36433468
O        6.30737167     1.91777319     0.31145472
C        4.47115922     3.41553138     0.60494183
H        5.09069398     4.28245626     0.75641911
C        3.12407502     3.49552153     0.63432307
H        2.60123483     4.42396853     0.80962128
N        2.32034427     2.40483955     0.44391704
H        1.29629244     2.47478724     0.46770730
C        2.82027675     1.15461676     0.20974482
O        2.10824430     0.16511187     0.03627464
units angstrom
}
""", 0)

S66_23 = input.process_input("""
molecule dimer {
0 1
C       -1.23272700     1.21163896    -0.14162406
O       -0.57127667     2.24201573     0.02561679
N       -0.67058051     0.00388878    -0.31428147
H        0.34384695    -0.09056011    -0.30832667
H       -1.24421373    -0.80632370    -0.44668271
C       -2.73824495     1.26675766    -0.15588657
H       -3.07797534     1.64660511     0.80450159
H       -3.20211503     0.30286549    -0.34621112
H       -3.04998747     1.97549049    -0.91859737
--
0 1
N        4.19521289     1.11742864    -0.11954193
H        4.68524234     0.24147146    -0.23748040
C        4.99883890     2.26027358     0.03093977
O        6.21440093     2.16465126     0.01575499
C        4.22624673     3.47559007     0.19408371
H        4.74800972     4.40878293     0.31711883
C        2.87708602     3.41391454     0.18840695
H        2.25668197     4.29027492     0.30608385
N        2.19200391     2.24163303     0.03384119
H        1.15921343     2.23257196     0.03300387
C        2.82289388     1.03716353    -0.12841885
O        2.22570515    -0.02675243    -0.27022634
units angstrom
}
""", 0)

S66_24 = input.process_input("""
molecule dimer {
0 1
C        0.71264532     1.12099570     0.06054078
H        1.35784165     1.98639917     0.12773717
C        1.25823573    -0.15925190     0.12423352
H        2.32495428    -0.28709988     0.24674303
C        0.42688496    -1.27452666     0.04265043
H        0.85044465    -2.26843268     0.09474995
C       -0.94957784    -1.11007406    -0.10031360
H       -1.59445570    -1.97627370    -0.16371348
C       -1.49552564     0.17105056    -0.16154602
H       -2.56378279     0.29922115    -0.27370311
C       -0.66382760     1.28664289    -0.08340143
H       -1.08690070     2.28100020    -0.13288613
--
0 1
C        1.98776046     1.10975720     3.71031958
H        2.63260558     1.97594094     3.77407030
C        2.53371358    -0.17139390     3.77183931
H        3.60192047    -0.29954095     3.88458353
C        1.70206410    -1.28699400     3.69318889
H        2.12514581    -2.28134643     3.74284255
C        0.32566254    -1.12135897     3.54847214
H       -0.31944006    -1.98676921     3.48083951
C       -0.21989733     0.15887378     3.48450631
H       -1.28652536     0.28670299     3.36132755
C        0.61137962     1.27415454     3.56657725
H        0.18785474     2.26805957     3.51420832
units angstrom
}
""", 0)

S66_25 = input.process_input("""
molecule dimer {
0 1
N        1.57248145     0.25454916    -0.25648131
C        0.96935990    -0.90316032     0.04452614
H        1.61363891    -1.77218120     0.10234520
C       -0.39815811    -1.02881911     0.28096043
H       -0.81842477    -1.99173710     0.53356364
C       -1.19580525     0.10655779     0.19539732
H       -2.26068964     0.04953865     0.37344280
C       -0.58712829     1.31741239    -0.12010544
H       -1.16181223     2.22950003    -0.20046257
C        0.78854733     1.33970567    -0.33224053
H        1.28843202     2.26879436    -0.57852690
--
0 1
N       -0.53372327    -1.51586163     3.84414371
C       -1.46620136    -0.55523217     3.91799487
H       -2.46899061    -0.88618697     4.16018773
C       -1.20419832     0.79583625     3.70861549
H       -2.00275608     1.52034169     3.78688658
C        0.09522901     1.18507754     3.39834708
H        0.33721357     2.22407602     3.22247582
C        1.07478832     0.20217938     3.31498561
H        2.09708956     0.44892512     3.06654863
C        0.71230860    -1.12295838     3.54817861
H        1.45616936    -1.90851301     3.49173001
units angstrom
}
""", 0)

S66_26 = input.process_input("""
molecule dimer {
0 1
N        1.37690111     0.83974747     0.73462494
H        1.05181240     1.38622385     1.52335563
C        1.30898271     1.45752981    -0.52065500
O        0.92056136     2.61107777    -0.62597673
N        2.01142293    -1.21320830    -0.09807182
H        1.72728551     0.99084268    -2.61199556
C        2.02573687    -0.69717123    -1.36439740
H        2.29751698    -1.39106004    -2.14564531
C        1.71451235     0.59193780    -1.61248722
H        2.12945422    -2.20152091     0.05682913
C        1.64594503    -0.48520598     1.01871830
O        1.56111602    -0.97181638     2.12980905
--
0 1
N       -1.35546089    -0.83604594     0.73462494
H       -1.03037218    -1.38252232     1.52335563
C       -1.28754249    -1.45382828    -0.52065500
O       -0.89912114    -2.60737623    -0.62597673
N       -1.98998271     1.21690983    -0.09807182
H       -1.70584529    -0.98714115    -2.61199556
C       -2.00429665     0.70087276    -1.36439740
H       -2.27607676     1.39476157    -2.14564531
C       -1.69307213    -0.58823627    -1.61248722
H       -2.10801399     2.20522244     0.05682913
C       -1.62450481     0.48890751     1.01871830
O       -1.53967580     0.97551791     2.12980905
units angstrom
}
""", 0)

S66_27 = input.process_input("""
molecule dimer {
0 1
C        0.81874699     0.86417234     0.18828612
H        1.46611361     1.71666767     0.34472141
C        1.36899712    -0.39052394    -0.06669818
H        2.44303637    -0.51186194    -0.11057444
C        0.53437860    -1.48849320    -0.27188804
H        0.96084825    -2.46156422    -0.47550749
C       -0.84911561    -1.33050735    -0.21989643
H       -1.49706942    -2.18186028    -0.37955321
C       -1.39948546    -0.07603020     0.04043417
H       -2.47268667     0.04490778     0.09338206
C       -0.56529230     1.02140336     0.24227921
H       -0.99255667     1.99366131     0.44625817
--
0 1
N       -2.39843199     0.16214088     3.52041137
C       -1.78354606     1.31980869     3.80047556
H       -2.43115011     2.17298014     3.96298765
C       -0.40133116     1.46065642     3.89064637
H        0.03051760     2.42430654     4.12186267
C        0.39962023     0.34367712     3.67643246
H        1.47718940     0.41406140     3.73126697
C       -0.22093167    -0.86497792     3.38277288
H        0.35484284    -1.76059980     3.19869795
C       -1.61144595    -0.90301580     3.31732347
H       -2.12029887    -1.83146918     3.08848079
units angstrom
}
""", 0)

S66_28 = input.process_input("""
molecule dimer {
0 1
C        0.82576911     1.23652484    -0.04025044
H        1.52101317     2.06312520    -0.08247145
C        1.30015992    -0.06294088     0.12725601
H        2.36365753    -0.24226113     0.20767420
C        0.40352312    -1.12855218     0.19824486
H        0.77375338    -2.13742677     0.32412109
C       -0.96780949    -0.89519049     0.10313994
H       -1.66520900    -1.71998342     0.16042745
C       -1.44350838     0.40448328    -0.06244130
H       -2.50751124     0.58550112    -0.12415016
C       -0.54575549     1.46876875    -0.13624741
H       -0.91422190     2.47742220    -0.26785516
--
0 1
N       -0.27488064     0.67158742     3.21864568
H       -0.64818803     1.57334885     2.95575271
C        1.11726604     0.59860052     3.35065902
O        1.80817636     1.59302421     3.20582496
C        1.59616616    -0.73547719     3.66876922
H        2.65321825    -0.88769313     3.80289036
C        0.71645693    -1.74985837     3.79498575
H        1.02238445    -2.75827898     4.03151011
N       -0.62878896    -1.56482645     3.62489361
H       -1.27753679    -2.32738539     3.72376278
C       -1.20323727    -0.34002542     3.32547899
O       -2.40102568    -0.18920215     3.18336680
units angstrom
}
""", 0)

S66_29 = input.process_input("""
molecule dimer {
0 1
N        1.21075533     0.02867578     0.32971111
C        0.61193497    -1.15844901     0.15345176
H        1.25147791    -2.02952340     0.21929295
C       -0.75131399    -1.30864956    -0.08883407
H       -1.17041577    -2.29686932    -0.21338320
C       -1.54786767    -0.16994027    -0.15646691
H       -2.61101275    -0.24595469    -0.33875574
C       -0.94362237     1.07063612     0.01982310
H       -1.51881431     1.98450028    -0.01164403
C        0.42771857     1.11610863     0.25734879
H        0.92469451     2.06805173     0.39754798
--
0 1
N       -0.71316758    -0.28394932     3.29752332
H       -1.60805660    -0.71581281     3.11291983
C       -0.71291270     1.11386048     3.39053432
O       -1.75279577     1.74206028     3.27568419
C        0.60658206     1.67294182     3.61809739
H        0.70789842     2.74016399     3.71396557
C        1.67645565     0.85424952     3.68961744
H        2.68033469     1.22291422     3.83804398
N        1.55839451    -0.50304375     3.57706278
H        2.37183050    -1.09523110     3.56889514
C        0.35794757    -1.15027617     3.35068108
O        0.26581032    -2.35569425     3.21710180
units angstrom
}
""", 0)

S66_30 = input.process_input("""
molecule dimer {
0 1
C        0.83551718     1.11516693     0.02140131
H        1.48432398     1.98060858     0.01953430
C        1.38327497    -0.16614721     0.02376531
H        2.45714902    -0.29520468     0.02277108
C        0.54755466    -1.28131632     0.02168563
H        0.97293610    -2.27580453     0.01977853
C       -0.83552313    -1.11516159     0.02139907
H       -1.48433419    -1.98060640     0.01953009
C       -1.38328358     0.16615413     0.02375775
H       -2.45715618     0.29520906     0.02275707
C       -0.54756577     1.28132347     0.02168025
H       -0.97294284     2.27580548     0.01976873
--
0 1
C        0.65578060    -0.11679048     3.53075174
H        1.04724138    -1.12390931     3.52628348
H        1.37085438     0.69327350     3.52625015
C       -0.65577592     0.11679215     3.53076063
H       -1.37084787    -0.69327237     3.52626454
H       -1.04723903     1.12391105     3.52630243
units angstrom
}
""", 0)

S66_31 = input.process_input("""
molecule dimer {
0 1
N       -0.05087365    -0.98008127     0.03396219
H       -0.05322205    -1.99069374     0.04982167
C       -1.30881316    -0.36187638     0.00402596
O       -2.32722000    -1.03255492    -0.00582886
C       -1.23681849     1.08804829    -0.01222440
H       -2.15273897     1.65146044    -0.05477443
C       -0.03519433     1.69783584     0.03370483
H        0.07036636     2.77247575     0.03188224
N        1.13452913     0.99028251     0.09184461
H        2.02372032     1.45677218     0.15569277
C        1.19318599    -0.39183287     0.11577512
O        2.23639797    -1.01118826     0.19418562
--
0 1
C        0.72600726     0.02505349     3.39819044
H        1.24312499    -0.84593440     3.02096384
H        1.33161826     0.81204754     3.82550477
C       -0.60276924     0.12564394     3.34894351
H       -1.21477213    -0.66183565     2.93204279
H       -1.11459423     0.99671353     3.73294327
units angstrom
}
""", 0)

S66_32 = input.process_input("""
molecule dimer {
0 1
N       -0.05545357    -0.94799090     0.01001028
H       -0.05731609    -1.95771330     0.05505287
C       -1.31395971    -0.33514498    -0.06458622
O       -2.32889664    -1.00790087    -0.12310273
C       -1.24835877     1.11605191    -0.06650860
H       -2.16434937     1.67533298    -0.14710244
C       -0.05308010     1.73142748     0.03419541
H        0.04811054     2.80642986     0.04341968
N        1.11592628     1.02759107     0.13516893
H        1.99665515     1.49727976     0.26162029
C        1.17534700    -0.35380470     0.17616616
O        2.21463146    -0.96646542     0.33517250
--
0 1
C        0.70785184    -0.17230221     3.27635136
H        1.70367011    -0.52628807     3.16213263
C       -0.43675225     0.21415547     3.38254320
H       -1.44163480     0.54285582     3.48290737
units angstrom
}
""", 0)

S66_33 = input.process_input("""
molecule dimer {
0 1
N        1.38138219    -0.00023348     0.13146374
C        0.67935079    -1.14023946     0.09207966
H        1.25871960    -2.05496223     0.12588361
C       -0.70972232    -1.19311407     0.00666426
H       -1.21408768    -2.14856163    -0.02530851
C       -1.42161357     0.00013343    -0.04081690
H       -2.50069615     0.00025757    -0.10916973
C       -0.70940120     1.19317538     0.00652198
H       -1.21351163     2.14874784    -0.02552831
C        0.67965167     1.13995623     0.09189303
H        1.25926073     2.05451090     0.12550248
--
0 1
C        0.01960458     0.66643934     3.48727228
H        0.93007858     1.22592506     3.32815744
H       -0.88994292     1.22884357     3.64423278
C        0.01993726    -0.66624796     3.48740452
H        0.93067296    -1.22533044     3.32839408
H       -0.88935083    -1.22907273     3.64449367
units angstrom
}
""", 0)

S66_34 = input.process_input("""
molecule dimer {
0 1
C       -2.53330865    -0.29487907     0.71314876
H       -2.56362682    -0.97708181    -0.13642264
H       -2.56697835    -0.89587590     1.62173177
H       -3.43442611     0.31595713     0.68410447
C       -1.27188487     0.55765547     0.67435468
H       -1.27102630     1.25656571     1.51431940
H       -1.26663255     1.16789581    -0.23182653
C       -0.00013504    -0.27841822     0.71960315
H       -0.00015938    -0.88722952     1.62863709
H       -0.00036543    -0.98071418    -0.11940439
C        1.27189476     0.55738219     0.67406108
H        1.27097175     1.25663331     1.51370541
H        1.26663649     1.16718250    -0.23238692
C        2.53340376    -0.29494176     0.71328015
H        2.56391919    -0.97777410    -0.13577836
H        3.43430956     0.31625432     0.68359945
H        2.56755821    -0.89520887     1.62232865
--
0 1
C        2.53355730     0.29502133     4.51309986
H        2.56814179     0.89482803     3.60377431
H        2.56406061     0.97822791     5.36184468
H        3.43423799    -0.31647598     4.54330880
C        1.27173110    -0.55686594     4.55240411
H        1.26628739    -1.16659365     5.45890107
H        1.27060059    -1.25621968     3.71282305
C       -0.00004389     0.27923316     4.50678767
H       -0.00019882     0.98154314     5.34577214
H        0.00003301     0.88800958     3.59771803
C       -1.27180473    -0.55690882     4.55205921
H       -1.26642249    -1.16701827     5.45830931
H       -1.27069839    -1.25593171     3.71219555
C       -2.53352396     0.29513749     4.51308150
H       -2.56771726     0.89567116     3.60420474
H       -3.43432593    -0.31616087     4.54259468
H       -2.56406349     0.97772373     5.36234289
units angstrom
}
""", 0)

S66_35 = input.process_input("""
molecule dimer {
0 1
C       -2.53038287    -0.41757533     0.68130643
H       -2.55988603    -0.98278998    -0.25015619
H       -2.55403625    -1.13386495     1.50265790
H       -3.43621355     0.18414376     0.73677133
C       -1.27615683     0.44363493     0.75002483
H       -1.27808384     1.02521785     1.67508548
H       -1.28033899     1.16855564    -0.06715806
C        0.00220470    -0.38071620     0.67899257
H        0.00782894    -1.11141304     1.49383122
H        0.00624866    -0.96052270    -0.24882046
C        1.26833347     0.46239635     0.74936913
H        1.26201986     1.04425029     1.67424645
H        1.26163488     1.18705711    -0.06803458
C        2.53496627    -0.38042469     0.68068636
H        2.57244024    -0.94571652    -0.25045186
H        3.43198117     0.23441492     0.73557772
H        2.56920771    -1.09581003     1.50245608
--
0 1
C       -0.00052120     0.06397129     5.24130633
C        0.00055054    -0.07615981     6.76103928
H       -0.88648549     0.38791623     7.19440870
H        0.00980204    -1.12694006     7.05404915
H        0.87921076     0.40350475     7.19468235
C       -1.23997654    -0.61768074     4.66740782
H       -1.26327576    -0.52872361     3.58057863
H       -1.25206217    -1.67895713     4.92042102
H       -2.15092026    -0.16538948     5.06249294
C        1.25208391    -0.59356951     4.66783599
H        1.27341069    -0.50528385     3.58086503
H        1.28521444    -1.65413035     4.92192831
H        2.15389614    -0.12292620     5.06225711
C       -0.01476908     1.54376378     4.86668505
H        0.86299692     2.05435080     5.26564018
H       -0.01529328     1.67021871     3.78303336
H       -0.90287503     2.03709750     5.26447319
units angstrom
}
""", 0)

S66_36 = input.process_input("""
molecule dimer {
0 1
C        0.38252221    -0.07060697     0.76689582
C       -1.04063947     0.39681125     1.06093593
H       -1.77157460    -0.28150025     0.61833023
H       -1.22471777     0.43573509     2.13551890
H       -1.21406603     1.39372444     0.65309065
C        0.59084747    -1.46681814     1.34797791
H        1.60291380    -1.82295000     1.15010285
H        0.43896858    -1.46674598     2.42828668
H       -0.10991906    -2.17868425     0.90931390
C        1.37826905     0.89843536     1.39914944
H        2.40439397     0.58544074     1.20073365
H        1.24378092     0.94597430     2.48070991
H        1.24837318     1.90502262     0.99895071
C        0.60196094    -0.11103419    -0.74309659
H        0.45921182     0.87703910    -1.18289819
H        1.61369399    -0.44345945    -0.97967210
H       -0.09953078    -0.79754982    -1.21922069
--
0 1
C       -0.37502842     0.06931363     5.96648833
C        1.04778403    -0.39965237     5.67308879
H        1.23222323    -0.43898152     4.59856833
H        1.77921818     0.27802046     6.11582437
H        1.22004770    -1.39665841     6.08120936
C       -0.58142523     1.46587516     5.38565786
H       -1.59338833     1.82286061     5.58250538
H        0.11949337     2.17694663     5.82537963
H       -0.42831602     1.46607177     4.30551550
C       -0.59532291     0.10948985     7.47634196
H       -1.60653907     0.44376683     7.71241515
H        0.10718954     0.79443888     7.95318018
H       -0.45475982    -0.87903049     7.91579370
C       -1.37149114    -0.89846403     5.33334194
H       -1.24256513    -1.90543941     5.73292091
H       -2.39738024    -0.58469117     5.53172979
H       -1.23678678    -0.94543842     4.25176527
units angstrom
}
""", 0)

S66_37 = input.process_input("""
molecule dimer {
0 1
C        0.79991408    -1.02205164     0.68773696
H        0.85355588    -1.12205101    -0.39801435
H        1.49140210    -1.74416936     1.11972040
C        1.11688700     0.42495279     1.09966205
H        1.83814230     0.89014504     0.43045256
H        1.55556959     0.43982464     2.09708356
C       -0.24455916     1.16568959     1.10297714
H       -0.25807760     2.00086313     0.40532333
H       -0.44880450     1.57699582     2.09098447
C       -1.29871418     0.10381191     0.73930899
H       -1.47356078     0.10524338    -0.33800545
H       -2.25673428     0.27804118     1.22715843
C       -0.64687993    -1.22006836     1.13630660
H       -1.12443918    -2.08762702     0.68299327
H       -0.68601864    -1.34528332     2.22022006
--
0 1
C        0.04984615     0.09420760     5.61627735
C       -0.04649805    -0.05787837     7.13191782
H        0.94604832    -0.07334458     7.58427505
H       -0.60542282     0.77000613     7.57035274
H       -0.55366275    -0.98654445     7.39726741
C        0.76389939     1.40111272     5.28065247
H        0.84541894     1.53461185     4.20097059
H        0.22042700     2.25580115     5.68615385
H        1.77150393     1.41176313     5.69888547
C       -1.35516567     0.11403225     5.01895782
H       -1.31823408     0.23122219     3.93510886
H       -1.93746520     0.94145581     5.42730374
H       -1.88506873    -0.81375459     5.24028712
C        0.83774596    -1.07927730     5.03893917
H        0.34252564    -2.02626804     5.25918232
H        0.93258913    -0.99209454     3.95580439
H        1.84246405    -1.11668194     5.46268763
units angstrom
}
""", 0)

S66_38 = input.process_input("""
molecule dimer {
0 1
C        0.95688019    -0.89184563     1.14195000
H        1.50456597    -1.27835762     0.28342019
H        1.42138447    -1.31477793     2.03102546
C        0.99094943     0.65850830     1.14550384
H        1.51059446     1.02309646     0.25994788
H        1.51625823     1.05981813     2.01053703
C       -0.47945194     1.10231879     1.10387910
H       -0.61626861     2.06487722     0.61356737
H       -0.87474223     1.18907144     2.11806960
C       -1.18210650    -0.05279656     0.39334575
H       -0.94888216    -0.02683030    -0.67380459
H       -2.26566452    -0.03356474     0.50127403
C       -0.53065958    -1.27488954     1.03930959
H       -0.69039061    -2.19702093     0.48299221
H       -0.95084939    -1.41541197     2.03674782
--
0 1
C       -1.13198517    -0.38391856     5.05596626
H       -1.46511966    -0.14721994     4.04338190
H       -1.93677357    -0.92701702     5.54895277
C        0.18162128    -1.17946347     5.00820507
H        0.23156623    -1.83720616     4.14207124
H        0.26190891    -1.81082110     5.89259036
C        1.31093651    -0.11675764     5.00880116
H        1.93220146    -0.17743649     4.11692754
H        1.96834600    -0.26664069     5.86420633
C        0.60076314     1.24491110     5.11666799
H        0.42089996     1.65340289     4.12066887
H        1.18114710     1.97931461     5.67264126
C       -0.74128932     0.91043867     5.76647985
H       -1.48095789     1.70295043     5.66159855
H       -0.60124939     0.71879862     6.83302881
units angstrom
}
""", 0)

S66_39 = input.process_input("""
molecule dimer {
0 1
C        0.76554546     0.86824433     0.82099095
H        1.43747647     1.68000664     1.06510281
C        1.23765260    -0.44283807     0.79388795
H        2.27575877    -0.64853808     1.01771141
C        0.37223723    -1.48853667     0.47726862
H        0.73818789    -2.50608012     0.45705609
C       -0.96493318    -1.22297162     0.18687834
H       -1.63645949    -2.03456079    -0.05777362
C       -1.43706509     0.08840558     0.21327714
H       -2.47468432     0.29430216    -0.01146746
C       -0.57190649     1.13402416     0.53081281
H       -0.93769935     2.15171058     0.55107764
--
0 1
C       -0.76345318    -0.72677383     4.05982770
H       -0.86970702    -0.55182467     2.98752083
H       -1.41509075    -1.55603772     4.33297836
C        0.70608801    -0.98383692     4.40395757
H        1.20131879    -1.62142197     3.67337330
H        0.76936719    -1.48405069     5.37142421
C        1.34622506     0.42155976     4.49491043
H        1.99649337     0.61423069     3.64305751
H        1.95909224     0.51072918     5.39063579
C        0.16717893     1.42073677     4.52178247
H        0.05002744     1.87970717     3.53949713
H        0.31277252     2.22224160     5.24418107
C       -1.06659283     0.56364158     4.81743133
H       -1.99758134     1.03937903     4.51151819
H       -1.13201859     0.35432067     5.88796657
units angstrom
}
""", 0)

S66_40 = input.process_input("""
molecule dimer {
0 1
C        0.31195353     0.56102334     0.49669886
H        0.74213608     1.55336911     0.48156571
C        1.14218235    -0.55807461     0.53606185
H        2.21651131    -0.43425014     0.55235015
C        0.58780415    -1.83668705     0.55414435
H        1.23191239    -2.70484153     0.58522179
C       -0.79665772    -1.99637562     0.53296300
H       -1.22677442    -2.98844427     0.54863708
C       -1.62689297    -0.87747365     0.49416828
H       -2.70112211    -1.00134997     0.47981498
C       -1.07266525     0.40120590     0.47597397
H       -1.71697357     1.26940117     0.44591995
--
0 1
C        0.17046797     0.50613197     4.83469402
C        1.61671665     0.68491933     4.37973254
H        2.03257337     1.61819721     4.76315552
H        2.24011597    -0.13569629     4.73858640
H        1.67732578     0.70431062     3.29079832
C        0.11607660     0.47476083     6.35955934
H       -0.90971343     0.34734041     6.70864711
H        0.71148250    -0.35092603     6.75211308
H        0.50437108     1.40264546     6.78246492
C       -0.37891207    -0.80336000     4.27439800
H       -1.41378567    -0.95363504     4.58706959
H        0.20754451    -1.65233376     4.63020927
H       -0.35013224    -0.80381278     3.18408376
C       -0.67090481     1.67070366     4.31848855
H       -0.64936386     1.70673405     3.22848999
H       -1.71069396     1.56693409     4.63297103
H       -0.29525222     2.62139813     4.70059546
units angstrom
}
""", 0)

S66_41 = input.process_input("""
molecule dimer {
0 1
N       -0.20890478    -0.96458262     0.53476104
H       -0.22415099    -1.97310940     0.60508386
C       -1.44634208    -0.34458112     0.30665858
O       -2.46123675    -1.01079161     0.19789196
C       -1.35778219     1.10318559     0.22814378
H       -2.25657214     1.66773071     0.04984731
C       -0.16300320     1.70989257     0.38112632
H       -0.04629046     2.78244591     0.33334968
N        0.98545210     1.00082412     0.61120636
H        1.86755978     1.46692777     0.74478430
C        1.02702092    -0.37917011     0.71264723
O        2.04919670    -0.99739548     0.93725979
--
0 1
C        1.14141247     2.35703152     4.05707817
H        0.71056385     2.66808022     3.10429560
H        0.50717856     2.76246464     4.84532582
H        2.12429249     2.81747894     4.15019966
C        1.21442893     0.83816057     4.14659651
H        1.64481257     0.54859772     5.10788747
H        1.88901852     0.44700002     3.38147835
C       -0.15035626     0.17999392     3.99177975
H       -0.82160052     0.54886973     4.77339899
H       -0.59782713     0.49025894     3.04187953
C       -0.09406732    -1.34069263     4.05141525
H        0.32953817    -1.64312304     5.01205144
H        0.59745442    -1.70257157     3.28691282
C       -1.46335024    -1.98256584     3.86764160
H       -1.90172924    -1.70910816     2.90745609
H       -1.40641145    -3.06933423     3.91169879
H       -2.15131302    -1.65421986     4.64687465
units angstrom
}
""", 0)

S66_42 = input.process_input("""
molecule dimer {
0 1
N        0.19572959    -0.84468925     0.82384642
H        0.45039753    -1.79675294     1.04976794
C       -1.17904919    -0.57368440     0.75948349
O       -1.99364624    -1.45626526     0.96690066
C       -1.47671471     0.81115567     0.43755952
H       -2.50635592     1.11565059     0.36389469
C       -0.46811280     1.68296245     0.23489084
H       -0.63843522     2.72164296    -0.00616410
N        0.84562854     1.30599113     0.32683051
H        1.58969256     1.96887924     0.18595979
C        1.25426147     0.01946187     0.63624397
O        2.42230438    -0.30171639     0.73187948
--
0 1
C        1.05672314    -0.86351031     4.39874366
H        1.51057565    -0.95556655     3.41076111
H        1.60122564    -1.52749058     5.06794134
C        1.11103661     0.60244169     4.83167965
H        2.06932660     1.07534062     4.62095536
H        0.92292133     0.68407923     5.90490278
C       -0.05631497     1.21525617     4.06090845
H        0.21798930     1.30403777     3.00743682
H       -0.34072939     2.20639729     4.41254246
C       -1.17325946     0.17768426     4.23193676
H       -1.89879874     0.20129811     3.42056485
H       -1.71734509     0.38238141     5.15418538
C       -0.45022312    -1.18886357     4.33559365
H       -0.69288766    -1.83301970     3.49223397
H       -0.76532935    -1.71626599     5.23468007
units angstrom
}
""", 0)

S66_43 = input.process_input("""
molecule dimer {
0 1
N        0.62608128    -0.85091265     0.80591569
H        0.40918989    -1.81150056     1.03440142
C       -0.43245619    -0.08733581     0.29466376
O       -1.53077162    -0.58840313     0.12359257
C       -0.06687462     1.29127521     0.01963739
H       -0.80974352     1.95181039    -0.39283965
C        1.18354208     1.71793501     0.29053321
H        1.50185022     2.73387064     0.10983284
N        2.13412979     0.88660160     0.81908177
H        3.05533594     1.22390137     1.04342778
C        1.90278319    -0.44317844     1.12831175
O        2.74380631    -1.16392354     1.62858730
--
0 1
C       -0.62370220    -0.02971796     4.73188916
C       -1.94044838     0.71157084     4.94676206
H       -2.64751979     0.09336465     5.50162440
H       -1.78094882     1.63175538     5.51094708
H       -2.39815816     0.97306786     3.99160840
C       -0.00826558    -0.38315588     6.08316660
H        0.93489659    -0.91552919     5.95238477
H        0.18875537     0.51658585     6.66796874
H       -0.67955960    -1.02089289     6.65990335
C        0.34142207     0.86375986     3.95610006
H        1.28999256     0.35116515     3.78574607
H        0.54671227     1.78189631     4.50952643
H       -0.08097331     1.14224647     2.98863562
C       -0.88501939    -1.30975236     3.94152426
H       -1.34875779    -1.08791865     2.97889962
H        0.04755691    -1.84815128     3.76188758
H       -1.55552720    -1.97156632     4.49170918
units angstrom
}
""", 0)

S66_44 = input.process_input("""
molecule dimer {
0 1
C        0.66640038     0.18381078     0.41973683
H        1.22888182    -0.32988301     1.18625971
H        1.22803556     0.69720813    -0.34760989
C       -0.66597358     0.18297343     0.41961191
H       -1.22792171    -0.33149890     1.18610334
H       -1.22818427     0.69564575    -0.34774808
--
0 1
C       -2.53275995    -0.39365922     4.14534248
H       -2.56225339    -1.00668000     3.24415261
H       -2.56889390    -1.06787984     5.00095950
H       -3.43393131     0.21735721     4.16258843
C       -1.27132347     0.45901620     4.18116042
H       -1.27172933     1.07910977     5.08055437
H       -1.26293512     1.14592451     3.33210001
C       -0.00004920    -0.37854138     4.15421721
H       -0.00020326    -1.06521408     5.00604923
H        0.00009186    -1.00611921     3.25757472
C        1.27117120     0.45904505     4.18162175
H        1.27144420     1.07885580     5.08110716
H        1.26297638     1.14611970     3.33271412
C        2.53262258    -0.39367946     4.14579757
H        2.56224605    -1.00653596     3.24448839
H        3.43380069     0.21725671     4.16337561
H        2.56854094    -1.06813554     5.00130328
units angstrom
}
""", 0)

S66_45 = input.process_input("""
molecule dimer {
0 1
C       -0.60618936     0.05587406     0.58900491
H       -1.66803667     0.05577624     0.58901162
C        0.60584873     0.05554087     0.58926624
H        1.66767817     0.05486328     0.58972794
--
0 1
C       -2.53040391    -0.34745600     4.21851416
H       -2.53877054    -1.00940954     3.35210357
H       -2.58232224    -0.97372522     5.10910493
H       -3.43281853     0.26144806     4.18575253
C       -1.26987178     0.50714472     4.22958343
H       -1.28652345     1.18014394     5.08999255
H       -1.24460479     1.14136072     3.34078732
C        0.00004684    -0.33118629     4.27003876
H        0.00004957    -0.94897593     5.17310016
H        0.00011393    -1.01948544     3.42079757
C        1.26994540     0.50718978     4.22967030
H        1.28657322     1.18015690     5.09009161
H        1.24480048     1.14136210     3.34086911
C        2.53046789    -0.34744680     4.21872389
H        2.53884766    -1.00942955     3.35234481
H        3.43284666     0.26148455     4.18599753
H        2.58228512    -0.97366153     5.10935743
units angstrom
}
""", 0)

S66_46 = input.process_input("""
molecule dimer {
0 1
C        1.37219093     1.01247736     0.97082468
H        0.95217623     2.01404955     1.03311725
H        1.94742170     0.92651560     0.05071776
H        2.05170208     0.85182517     1.80295247
C        0.32673706    -0.07764727     0.98819876
O        0.61882128    -1.25248130     1.17128126
N       -0.95002884     0.34488680     0.77391491
H       -1.10467156     1.32202550     0.60611216
C       -2.05985440    -0.57736895     0.68015349
H       -1.66935602    -1.56679601     0.89718425
H       -2.83459176    -0.33138032     1.40366139
H       -2.49097050    -0.57892483    -0.31993926
--
0 1
C        2.66066552     0.46274539     4.85334645
H        2.77750480     1.21716129     4.07460163
H        2.57455515     0.98763172     5.80500251
H        3.57275696    -0.13149652     4.88015446
C        1.43239329    -0.40064212     4.59579490
H        1.33782394    -1.14609612     5.38884574
H        1.54881342    -0.95410645     3.66195110
C        0.14985545     0.41797183     4.53049355
H        0.03828513     0.99570671     5.45357719
H        0.22908959     1.15078674     3.72084090
C       -1.09450084    -0.43236340     4.31361365
H       -1.18530281    -1.14684989     5.13503088
H       -0.96669384    -1.02130113     3.40339920
C       -2.36133934     0.40792810     4.22349893
H       -2.29442610     1.11497908     3.39572969
H       -3.24668156    -0.20808939     4.06966602
H       -2.51169538     0.98413919     5.13671852
units angstrom
}
""", 0)

S66_47 = input.process_input("""
molecule dimer {
0 1
C        0.72918867     1.11310122     0.32672825
H        1.30321590     2.01422234     0.15916027
C        1.37508737    -0.11936635     0.41277695
H        2.45051474    -0.17462400     0.31330720
C        0.63503981    -1.28055339     0.62938541
H        1.13633448    -2.23601747     0.70021716
C       -0.75098563    -1.20965430     0.75789034
H       -1.32452590    -2.11141283     0.92419891
C       -1.39703443     0.02267081     0.67308963
H       -2.47242537     0.07848826     0.77399799
C       -0.65689731     1.18429622     0.45833859
H       -1.15782845     2.14058713     0.39509608
--
0 1
C        0.15810619     0.15289032     4.08588285
H        0.28023260     0.37837378     3.03545641
C       -0.93297537    -0.60200829     4.51321912
H       -1.65347990    -0.95852255     3.78952470
C       -1.09367536    -0.89613361     5.86616918
H       -1.94078294    -1.48210218     6.19641672
C       -0.16179279    -0.43508023     6.79466326
H       -0.28568629    -0.66304639     7.84467076
C        0.92979230     0.32002182     6.36942298
H        1.65291139     0.67785500     7.08980563
C        1.08859620     0.61350684     5.01593166
H        1.93585412     1.19958163     4.68588434
units angstrom
}
""", 0)

S66_48 = input.process_input("""
molecule dimer {
0 1
N        1.32276272    -0.01037598     1.01918373
C        0.65128601    -1.14899203     0.79680119
H        1.20041842    -2.06552808     0.97367282
C       -0.67268130    -1.19471172     0.36665693
H       -1.15719362    -2.14732141     0.20646407
C       -1.34719676     0.00313399     0.15214401
H       -2.37535653     0.00840542    -0.18229302
C       -0.66455797     1.19409062     0.37900199
H       -1.14262633     2.15155765     0.22872051
C        0.65889576     1.13497854     0.80885987
H        1.21410272     2.04591045     0.99543831
--
0 1
N        0.45011507     0.00130104     6.78095972
C        1.32078309    -0.00431175     5.76154669
H        2.36863966    -0.00306323     6.03584948
C        0.94739735    -0.01137951     4.41971862
H        1.69485802    -0.01554353     3.63861897
C       -0.40865120    -0.01279358     4.10730315
H       -0.73837988    -0.01824905     3.07702170
C       -1.32675447    -0.00707849     5.15247277
H       -2.39120450    -0.00792788     4.96373698
C       -0.85115066    -0.00016084     6.46143162
H       -1.54333433     0.00442229     7.29462282
units angstrom
}
""", 0)

S66_49 = input.process_input("""
molecule dimer {
0 1
C        0.84507720     1.05791869     0.69945490
H        1.50640601     1.90322178     0.83338235
C        1.37550931    -0.21745534     0.51116093
H        2.44718367    -0.36147258     0.50285232
C        0.52406810    -1.30704432     0.33319233
H        0.93572726    -2.29602641     0.18492305
C       -0.85771573    -1.12146341     0.34638409
H       -1.51838119    -1.96645805     0.20836325
C       -1.38804570     0.15363438     0.53761349
H       -2.45971752     0.29741587     0.55003229
C       -0.53661315     1.24342221     0.71273882
H       -0.94892427     2.23280628     0.85736635
--
0 1
N        0.02311730     0.35202455     6.77454464
C        0.17780112     1.28998616     5.82966776
H        0.31957195     2.30251216     6.18756949
C        0.16359185     1.02269639     4.46316833
H        0.29383191     1.82372219     3.74928292
C       -0.02074646    -0.28893329     4.03787790
H       -0.03731291    -0.53205196     2.98452996
C       -0.18259538    -1.27396762     5.00673698
H       -0.32913840    -2.30917859     4.73196547
C       -0.15339291    -0.90663452     6.34982649
H       -0.27698904    -1.65414849     7.12392749
units angstrom
}
""", 0)

S66_50 = input.process_input("""
molecule dimer {
0 1
C        0.83661195     1.11485600     0.23100790
H        1.48545250     1.97968049     0.21470491
C        1.38418781    -0.16696533     0.26005688
H        2.45768419    -0.29628753     0.26605977
C        0.54747934    -1.28184652     0.28693051
H        0.97191784    -2.27597918     0.31387670
C       -0.83666710    -1.11500365     0.28456279
H       -1.48555353    -1.97956851     0.30969784
C       -1.38416274     0.16685015     0.25560540
H       -2.45764469     0.29645927     0.25854055
C       -0.54749833     1.28174826     0.22897743
H       -0.97214124     2.27600137     0.21116093
--
0 1
C        0.00585466     0.07515017     3.77945155
H        0.00284553     0.05759463     2.71537604
C        0.00951511     0.09473103     4.99182772
H        0.01262752     0.11190396     6.05302473
units angstrom
}
""", 0)

S66_51 = input.process_input("""
molecule dimer {
0 1
C       -0.60172996    -0.02857012     0.38493492
H       -1.66373543    -0.02852657     0.37901431
C        0.61010917    -0.02866364     0.38816379
H        1.67213544    -0.02879308     0.38796752
--
0 1
C       -0.00735998     0.10033739     4.14281190
H       -0.00396560     0.06660234     3.07951502
C       -0.01129640     0.13862741     5.35427728
H       -0.01456263     0.17200329     6.41518870
units angstrom
}
""", 0)

S66_52 = input.process_input("""
molecule dimer {
0 1
C        0.96408039     0.87509331     0.37801364
H        1.65982961     1.69993082     0.44604227
C        1.43105709    -0.41313344     0.11899152
H        2.48952453    -0.58720917    -0.01701261
C        0.53412766    -1.47763890     0.04241755
H        0.89696129    -2.47738839    -0.15201199
C       -0.83032682    -1.25360409     0.22085611
H       -1.52576001    -2.07962435     0.16411655
C       -1.29758715     0.03441261     0.48024263
H       -2.35439607     0.20801612     0.62856096
C       -0.40044509     1.09977921     0.56160137
H       -0.76045514     2.09376880     0.78475698
--
0 1
C       -0.11985517     0.53438939     4.36008118
O       -0.58804476     1.58383601     3.98082079
O        0.28335741    -0.44317387     3.52079591
H        0.11465259    -0.11726029     2.61939066
C        0.09009913     0.13740231     5.79148697
H       -0.21986702     0.94673889     6.44147585
H       -0.48598160    -0.75922167     6.00843808
H        1.13859655    -0.09872978     5.95650555
units angstrom
}
""", 0)

S66_53 = input.process_input("""
molecule dimer {
0 1
C        0.85556074     0.35853244     1.04975426
H        1.51382550     0.90267956     1.71276582
C        1.34289713    -0.67537866     0.25115740
H        2.39288384    -0.93334472     0.28196305
C        0.47780661    -1.37670110    -0.58781577
H        0.85608399    -2.17890753    -1.20682428
C       -0.87482983    -1.04255615    -0.63045178
H       -1.54540573    -1.58570014    -1.28241614
C       -1.36239729    -0.00701391     0.16584645
H       -2.41157102     0.25346723     0.13077885
C       -0.49844404     0.69315695     1.00699199
H       -0.86611090     1.49033989     1.63803696
--
0 1
C        0.08192937     0.49753072     4.80472861
O        0.32841872     1.54095697     4.21748933
N       -0.22211788    -0.65747581     4.15356127
H       -0.19691756    -0.66449114     3.14692466
H       -0.37789436    -1.51296813     4.64926298
C        0.10477407     0.40263889     6.31314609
H        1.13648787     0.48685118     6.64821988
H       -0.31712984    -0.52400410     6.69417176
H       -0.44469059     1.24648520     6.71991660
units angstrom
}
""", 0)

S66_54 = input.process_input("""
molecule dimer {
0 1
C        0.78014717    -0.60991473    -1.20755689
H        0.89619160    -1.13763959    -2.14414463
C        0.47794275     0.75099363    -1.20789541
H        0.35696423     1.27816780    -2.14405407
C        0.32728928     1.43186787    -0.00000000
H        0.09146503     2.48713922     0.00000000
C        0.47794275     0.75099363     1.20789541
H        0.35696423     1.27816780     2.14405407
C        0.78014717    -0.60991473     1.20755689
H        0.89619160    -1.13763959     2.14414463
C        0.93164831    -1.28998134     0.00000000
H        1.16848573    -2.34521369    -0.00000000
--
0 1
O       -2.74383121    -0.26926257     0.00000000
H       -2.57902721    -1.21398410     0.00000000
H       -1.85653027     0.10232776     0.00000000
units angstrom
}
""", 0)

S66_55 = input.process_input("""
molecule dimer {
0 1
C        0.75974918     1.03127506     0.37377239
H        1.43501626     1.87566427     0.37470462
C        1.26661779    -0.26736234     0.42127308
H        2.33491597    -0.42918019     0.45943234
C        0.39532054    -1.35599116     0.42490511
H        0.78866193    -2.36249259     0.46303549
C       -0.98220564    -1.14665441     0.38127024
H       -1.65765632    -1.99114019     0.38512100
C       -1.48934612     0.15114979     0.33757234
H       -2.55794704     0.31375049     0.30771900
C       -0.61877516     1.24033121     0.33388373
H       -1.01176161     2.24710690     0.30436922
--
0 1
O        0.04701895     0.30618537     3.68511328
H        0.13311917     0.35605847     2.72791973
C       -0.84913165    -0.75142870     3.96816832
H       -0.94485234    -0.80816328     5.04910445
H       -1.84128123    -0.57973096     3.54437811
H       -0.48267133    -1.71446977     3.60525680
units angstrom
}
""", 0)

S66_56 = input.process_input("""
molecule dimer {
0 1
C        0.69231523     1.08829204     0.32484124
H        1.28194880     1.99194678     0.25251578
C        1.31818722    -0.15687008     0.28689607
H        2.39314337    -0.21947636     0.18840681
C        0.55801841    -1.32195045     0.38139986
H        1.04391922    -2.28757380     0.35761542
C       -0.82755236    -1.24142187     0.51168501
H       -1.41670095    -2.14525152     0.58533927
C       -1.45341138     0.00367145     0.54838107
H       -2.52823255     0.06570272     0.64984254
C       -0.69346094     1.16840108     0.45622907
H       -1.17873534     2.13440989     0.48572685
--
0 1
N        0.27506479    -0.22271725     3.85890709
H        0.40968315    -0.17867675     2.85583573
H        0.41655736     0.72242949     4.19137936
C       -1.10103469    -0.62910066     4.13634288
H       -1.25891125    -0.65764767     5.21289841
H       -1.87233687     0.01128013     3.69622388
H       -1.25572667    -1.63866846     3.76072118
units angstrom
}
""", 0)

S66_57 = input.process_input("""
molecule dimer {
0 1
C        0.40877989     1.05102502     0.37553605
H        1.01193875     1.94854570     0.36807788
C        1.01916788    -0.19976963     0.28905343
H        2.09557130    -0.27183333     0.21719099
C        0.24172263    -1.35688270     0.29668995
H        0.71521633    -2.32658869     0.22807218
C       -1.14617971    -1.26425757     0.39390198
H       -1.74918186    -2.16192663     0.39940980
C       -1.75727780    -0.01396023     0.48295173
H       -2.83351378     0.05824368     0.55903918
C       -0.97968602     1.14420653     0.47228370
H       -1.45405142     2.11400088     0.53713589
--
0 1
C        0.24562178     1.95675759     4.25663541
H       -0.11252332     2.12248844     3.24334264
H        1.27020534     2.31346716     4.33807692
H       -0.35847510     2.53039342     4.95498813
C        0.20877544     0.50359448     4.67234424
O        0.49340385     0.15123306     5.81088230
N       -0.16361983    -0.36212226     3.69310315
H       -0.32474773    -0.00413152     2.76703481
C       -0.20041270    -1.78900149     3.91119021
H       -0.12232513    -1.95590903     4.98118644
H       -1.13565324    -2.20735207     3.54445210
H        0.62871378    -2.29287426     3.41385278
units angstrom
}
""", 0)

S66_58 = input.process_input("""
molecule dimer {
0 1
N       -0.94121124     0.79004136     0.01171891
C       -0.92275524    -0.55237814     0.03537875
H        0.05724051    -1.01558800     0.05135491
C       -2.07651907    -1.33301813     0.03929035
H       -1.99652895    -2.41058573     0.05887720
C       -3.31631294    -0.70333955     0.01759905
H       -4.23157489    -1.27908429     0.01979377
C       -3.34889528     0.68701881    -0.00708596
H       -4.28544414     1.22610455    -0.02465899
C       -2.14310382     1.38263356    -0.00889005
H       -2.13809974     2.46565258    -0.02778297
--
0 1
N        2.53321129    -0.95002930     0.04251789
C        3.73499010    -1.54320554     0.04459773
H        3.72976625    -2.62616799     0.06648690
C        4.94092634    -0.84824698     0.02059635
H        5.87736466    -1.38778216     0.02369036
C        4.90860873     0.54205748    -0.00715036
H        5.82398367     1.11730853    -0.02633187
C        3.66892840     1.17234361    -0.00962746
H        3.58915567     2.24990219    -0.03071603
C        2.51501483     0.39233399     0.01556620
H        1.53510443     0.85599657     0.01390336
units angstrom
}
""", 0)

S66_59 = input.process_input("""
molecule dimer {
0 1
C       -1.00686722    -0.03056821    -0.02477285
H        0.05900333    -0.06093974    -0.04936562
C       -2.21874380     0.00317347     0.00259920
H       -3.27927730     0.03352491     0.02720048
--
0 1
O        2.26390460    -0.14557006    -0.11547082
H        2.83426102    -0.73533944     0.38155611
H        2.83590044     0.20541797    -0.80084297
units angstrom
}
""", 0)

S66_60 = input.process_input("""
molecule dimer {
0 1
C       -0.61056257     0.22750310    -0.17060207
H        0.10738506     0.86143603    -0.63420924
C       -1.38627573    -0.52532550     0.37997353
H       -2.08070324    -1.17406739     0.85437937
--
0 1
C        2.83444960    -0.64143137     0.46593603
O        2.58027054     0.31467087    -0.23290172
O        1.88654498    -1.41577160     1.03362263
H        1.02554559    -1.04847261     0.76585149
C        4.21008475    -1.12288120     0.81608694
H        4.94847057    -0.48533112     0.34523661
H        4.33629527    -1.11102648     1.89612226
H        4.33236190    -2.15072575     0.48285261
units angstrom
}
""", 0)

S66_61 = input.process_input("""
molecule dimer {
0 1
C       -2.27534498    -0.13507494     0.83133387
H       -2.49071776    -0.72792669    -0.05756635
H       -2.22632382    -0.81844641     1.67882341
H       -3.11202566     0.54494342     0.98740008
C       -0.96169812     0.61927789     0.66939920
H       -0.78869920     1.25043181     1.54470266
H       -1.02617687     1.29544524    -0.18645838
C        0.22650217    -0.31471031     0.47998579
H        0.30944439    -0.97513911     1.34803794
H        0.03915056    -0.96599875    -0.37878983
C        1.54300168     0.42117452     0.26899951
H        1.71163863     1.10777177     1.10244654
H        1.46609466     1.04374331    -0.62529358
C        2.72757633    -0.52686091     0.13745931
H        2.58874155    -1.20321391    -0.70575734
H        3.66150100     0.01169308    -0.01596863
H        2.83519407    -1.13740994     1.03407512
--
0 1
C       -0.48356149    -0.28786315     4.12125154
O       -0.90617543    -1.40304340     3.92410496
O       -1.29725385     0.77110237     4.35384102
H       -2.19801596     0.41672183     4.31330528
C        0.95670557     0.12180293     4.13845692
H        1.58252864    -0.74837801     3.98030176
H        1.13274299     0.85607656     3.35533234
H        1.19401682     0.59110388     5.09025931
units angstrom
}
""", 0)

S66_62 = input.process_input("""
molecule dimer {
0 1
C       -2.58777605    -0.32310566     0.46945828
H       -2.61038910    -0.87636604    -0.46961946
H       -2.65974410    -1.05188654     1.27771411
H       -3.47603507     0.30562460     0.50896129
C       -1.30955982     0.49739424     0.58506260
H       -1.31725060     1.08326190     1.50634108
H       -1.26237673     1.21557375    -0.23677617
C       -0.05682966    -0.36826029     0.55844017
H       -0.08617526    -1.07335882     1.39587537
H       -0.05380919    -0.97684333    -0.35147393
C        1.23159606     0.44006559     0.63203246
H        1.21328340     1.05356193     1.53459305
H        1.26629733     1.13137662    -0.21310563
C        2.47257523    -0.44314441     0.61922148
H        2.52071888    -1.03526342    -0.29489695
H        3.38773437     0.14408974     0.68390871
H        2.45929703    -1.13936423     1.45861821
--
0 1
C        0.04216222     0.20124208     4.11650819
O        0.06907449     1.38631556     3.82466701
N        1.17474249    -0.55063556     4.21932814
H        2.04568275    -0.12805505     3.95066588
H        1.13580453    -1.54252223     4.35075106
C       -1.24805876    -0.53769541     4.38096202
H       -1.10080876    -1.49841677     4.86808639
H       -1.75428629    -0.69600434     3.43014867
H       -1.88600271     0.08954102     4.99623387
units angstrom
}
""", 0)

S66_63 = input.process_input("""
molecule dimer {
0 1
C        0.60678496     1.33042185     0.31643451
H        1.24649846     2.20226434     0.33035231
C        1.11808466     0.08724886     0.68511652
H        2.15005753    -0.00388678     0.99375824
C        0.29290229    -1.03608737     0.66910727
H        0.68849686    -2.00096149     0.95537797
C       -1.04283174    -0.91671112     0.28818964
H       -1.68270956    -1.78848825     0.27934903
C       -1.55358838     0.32734899    -0.07994317
H       -2.58923495     0.42028908    -0.37734619
C       -0.72804164     1.45084316    -0.06684834
H       -1.12362379     2.41565865    -0.35386143
--
0 1
C        0.41898688    -0.27167884     4.02497697
O        1.61447955    -0.10772809     4.10149274
O       -0.16051479    -1.48308380     4.22441532
H        0.57393607    -2.08419229     4.41745344
C       -0.60289735     0.77225268     3.70429579
H       -0.12460293     1.74319903     3.65747301
H       -1.05569745     0.53905649     2.74158774
H       -1.38774836     0.76671618     4.45679527
units angstrom
}
""", 0)

S66_64 = input.process_input("""
molecule dimer {
0 1
C        1.62971482     0.50301252     0.27011189
H        1.64157338     1.45923792    -0.24808286
H        2.31531919    -0.18355470    -0.21758635
H        1.96974564     0.64936024     1.29398105
C        0.26182776    -0.13286122     0.31456221
O        0.09925265    -1.30961602     0.61183995
N       -0.77350225     0.70251214     0.02207590
H       -0.56901138     1.66655677    -0.16581434
C       -2.15001214     0.26596865     0.09505328
H       -2.14473761    -0.81940745     0.10091210
H       -2.64054318     0.61582035     1.00360442
H       -2.70774393     0.62075110    -0.76826057
--
0 1
C       -0.04575608     0.51799706     3.77621664
H       -0.05063764     1.26017087     4.56209922
H       -0.69428883     0.68576570     2.92753308
C        0.72275422    -0.56896486     3.84602626
H        1.36805919    -0.74079051     4.69615412
H        0.71764224    -1.30416499     3.05371698
units angstrom
}
""", 0)

S66_65 = input.process_input("""
molecule dimer {
0 1
N       -0.08303249     0.00071459     1.05519999
C       -0.20285376    -1.14172585     0.36493369
H       -0.09848563    -2.05509795     0.93743262
C       -0.44678144    -1.19176367    -1.00451226
H       -0.53364921    -2.14585511    -1.50417155
C       -0.57468209     0.00343953    -1.70430948
H       -0.76368391     0.00448010    -2.76872670
C       -0.45345675     1.19724254    -1.00091647
H       -0.54563080     2.15227264    -1.49779508
C       -0.20931111     1.14450759     0.36836730
H       -0.11016707     2.05669726     0.94357396
--
0 1
C        0.47183602    -0.00605819     5.54171896
H        0.58724607    -0.00548400     6.59673278
C        0.33976626    -0.00660792     4.33547166
H        0.22161814    -0.00634549     3.27096619
units angstrom
}
""", 0)

S66_66 = input.process_input("""
molecule dimer {
0 1
N       -0.54105920     0.02957620    -0.20899508
H        0.05555335    -0.78611810    -0.13029335
H       -1.46966940    -0.27470845     0.05314338
C       -0.07879927     1.04239036     0.73845886
H       -0.72015294     1.91941377     0.67198026
H       -0.05075819     0.72382293     1.78551453
H        0.92643072     1.35660379     0.46199919
--
0 1
N        2.34185022    -1.25680010     0.03015300
C        2.68028654    -0.44445604    -0.98155948
H        2.13761932    -0.58899402    -1.90694084
C        3.65161580     0.54767776    -0.88119247
H        3.87646824     1.17201804    -1.73404317
C        4.31245587     0.71721920     0.33107196
H        5.07030981     1.47945653     0.44745609
C        3.97232296    -0.11774333     1.39019492
H        4.45491136    -0.02728109     2.35289557
C        2.98854139    -1.08253234     1.19101154
H        2.70245706    -1.74627994     1.99762219
units angstrom
}
""", 0)

# <<< Geometry Specification Strings >>>
rxnpattern = re.compile(r'^(.+)-(.+)-(.+)$')
GEOS = {}
for rxn in HRXN:

    GEOS['%s-%s-dimer'      % (dbse, rxn)] = eval('%s_%s' % (dbse, rxn))
    GEOS['%s-%s-monoA-CP'   % (dbse, rxn)] = str(eval('%s_%s' % (dbse, rxn))) + monoA_CP
    GEOS['%s-%s-monoB-CP'   % (dbse, rxn)] = str(eval('%s_%s' % (dbse, rxn))) + monoB_CP
    GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = str(eval('%s_%s' % (dbse, rxn))) + monoA_unCP
    GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = str(eval('%s_%s' % (dbse, rxn))) + monoB_unCP
