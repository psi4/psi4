"""
| Database of <description of members and reference energy type>.
| Geometries from <Reference>.
| Reference interaction energies from <Reference>.

- **benchmark**

  - ``'<benchmark_name>'`` <Reference>.
  - |dl| ``'<default_benchmark_name>'`` |dr| <Reference>.

- **subset**

  - ``'small'``
  - ``'large'``
  - ``'<subset>'`` <members_description>

"""
import re
import input

# <<< BENCH12 Database Module >>>
dbse = 'BENCH12'

# <<< Database Members >>>
HRXN = ['acene1', 'acene2', 'alka10', 'alka2', 'alka4', 'alka6', 'alka8', 'alke10', 'alke2', 'alke4', 'alke6', 'alke8', 'Ar', 'bag', 'boat1', 'boat2', 'book1', 'book2', 'cage', 'cyclic', 'h2o', 'h2o2', 'h2o3', 'h2o4', 'h2o5', 'h2o6', 'Kr', 'Ne', 'prism', 'S22_DD_bz_2_pd', 'S22_DD_bz_2_t', 'S22_DD_bz_me', 'S22_DD_ch4_2', 'S22_HB_formamide_2', 'S22_HB_formic_2', 'S22_HB_h2o_2', 'S22_HB_nh3_2', 'S22_MX_bz_h2o', 'S22_MX_bz_hcn', 'S22_MX_bz_nh3', 'S22_MX_c2h2_c2h4', 'thio1', 'thio2', ]
HRXN_SM = ['h2o', 'Ne']
HRXN_LG = ['alka10']
HRXN_alkenes = ['alke2','alke4','alke6','alke8','alke10']
HRXN_alkanes = ['alka2','alka4','alka6','alka8','alka10']
HRXN_acenes = ['acene1','acene2']
HRXN_h2o_size = ['h2o','h2o2','h2o3','h2o4','h2o5','h2o6']
HRXN_h2o_shape = ['bag', 'boat1', 'boat2', 'book1', 'book2', 'cage', 'cyclic', 'prism']
HRXN_S22_MX = ['S22_MX_bz_h2o', 'S22_MX_bz_hcn', 'S22_MX_bz_nh3', 'S22_MX_c2h2_c2h4']
HRXN_S22_DD = ['S22_DD_bz_2_pd', 'S22_DD_bz_2_t', 'S22_DD_bz_me', 'S22_DD_ch4_2']
HRXN_S22_HB = ['S22_HB_formamide_2', 'S22_HB_formic_2', 'S22_HB_h2o_2', 'S22_HB_nh3_2']
HRXN_atoms = ['Ne', 'Ar', 'Kr']

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV['%s-%s'            % (dbse, 'acene1'                )] = ['%s-%s-reagent'      % (dbse, 'acene1')]
RXNM['%s-%s'            % (dbse, 'acene1'                )] = dict(zip(ACTV['%s-%s' % (dbse, 'acene1')], [+1]))

ACTV['%s-%s'            % (dbse, 'acene2'                )] = ['%s-%s-reagent'      % (dbse, 'acene2')]
RXNM['%s-%s'            % (dbse, 'acene2'                )] = dict(zip(ACTV['%s-%s' % (dbse, 'acene2')], [+1]))

ACTV['%s-%s'            % (dbse, 'alka10'                )] = ['%s-%s-reagent'      % (dbse, 'alka10')]
RXNM['%s-%s'            % (dbse, 'alka10'                )] = dict(zip(ACTV['%s-%s' % (dbse, 'alka10')], [+1]))

ACTV['%s-%s'            % (dbse, 'alka2'                 )] = ['%s-%s-reagent'      % (dbse, 'alka2')]
RXNM['%s-%s'            % (dbse, 'alka2'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'alka2')], [+1]))

ACTV['%s-%s'            % (dbse, 'alka4'                 )] = ['%s-%s-reagent'      % (dbse, 'alka4')]
RXNM['%s-%s'            % (dbse, 'alka4'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'alka4')], [+1]))

ACTV['%s-%s'            % (dbse, 'alka6'                 )] = ['%s-%s-reagent'      % (dbse, 'alka6')]
RXNM['%s-%s'            % (dbse, 'alka6'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'alka6')], [+1]))

ACTV['%s-%s'            % (dbse, 'alka8'                 )] = ['%s-%s-reagent'      % (dbse, 'alka8')]
RXNM['%s-%s'            % (dbse, 'alka8'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'alka8')], [+1]))

ACTV['%s-%s'            % (dbse, 'alke10'                )] = ['%s-%s-reagent'      % (dbse, 'alke10')]
RXNM['%s-%s'            % (dbse, 'alke10'                )] = dict(zip(ACTV['%s-%s' % (dbse, 'alke10')], [+1]))

ACTV['%s-%s'            % (dbse, 'alke2'                 )] = ['%s-%s-reagent'      % (dbse, 'alke2')]
RXNM['%s-%s'            % (dbse, 'alke2'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'alke2')], [+1]))

ACTV['%s-%s'            % (dbse, 'alke4'                 )] = ['%s-%s-reagent'      % (dbse, 'alke4')]
RXNM['%s-%s'            % (dbse, 'alke4'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'alke4')], [+1]))

ACTV['%s-%s'            % (dbse, 'alke6'                 )] = ['%s-%s-reagent'      % (dbse, 'alke6')]
RXNM['%s-%s'            % (dbse, 'alke6'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'alke6')], [+1]))

ACTV['%s-%s'            % (dbse, 'alke8'                 )] = ['%s-%s-reagent'      % (dbse, 'alke8')]
RXNM['%s-%s'            % (dbse, 'alke8'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'alke8')], [+1]))

ACTV['%s-%s'            % (dbse, 'Ar'                    )] = ['%s-%s-reagent'      % (dbse, 'Ar')]
RXNM['%s-%s'            % (dbse, 'Ar'                    )] = dict(zip(ACTV['%s-%s' % (dbse, 'Ar')], [+1]))

ACTV['%s-%s'            % (dbse, 'bag'                   )] = ['%s-%s-reagent'      % (dbse, 'bag')]
RXNM['%s-%s'            % (dbse, 'bag'                   )] = dict(zip(ACTV['%s-%s' % (dbse, 'bag')], [+1]))

ACTV['%s-%s'            % (dbse, 'boat1'                 )] = ['%s-%s-reagent'      % (dbse, 'boat1')]
RXNM['%s-%s'            % (dbse, 'boat1'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'boat1')], [+1]))

ACTV['%s-%s'            % (dbse, 'boat2'                 )] = ['%s-%s-reagent'      % (dbse, 'boat2')]
RXNM['%s-%s'            % (dbse, 'boat2'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'boat2')], [+1]))

ACTV['%s-%s'            % (dbse, 'book1'                 )] = ['%s-%s-reagent'      % (dbse, 'book1')]
RXNM['%s-%s'            % (dbse, 'book1'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'book1')], [+1]))

ACTV['%s-%s'            % (dbse, 'book2'                 )] = ['%s-%s-reagent'      % (dbse, 'book2')]
RXNM['%s-%s'            % (dbse, 'book2'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'book2')], [+1]))

ACTV['%s-%s'            % (dbse, 'cage'                  )] = ['%s-%s-reagent'      % (dbse, 'cage')]
RXNM['%s-%s'            % (dbse, 'cage'                  )] = dict(zip(ACTV['%s-%s' % (dbse, 'cage')], [+1]))

ACTV['%s-%s'            % (dbse, 'cyclic'                )] = ['%s-%s-reagent'      % (dbse, 'cyclic')]
RXNM['%s-%s'            % (dbse, 'cyclic'                )] = dict(zip(ACTV['%s-%s' % (dbse, 'cyclic')], [+1]))

ACTV['%s-%s'            % (dbse, 'h2o'                   )] = ['%s-%s-reagent'      % (dbse, 'h2o')]
RXNM['%s-%s'            % (dbse, 'h2o'                   )] = dict(zip(ACTV['%s-%s' % (dbse, 'h2o')], [+1]))

ACTV['%s-%s'            % (dbse, 'h2o2'                  )] = ['%s-%s-reagent'      % (dbse, 'h2o2')]
RXNM['%s-%s'            % (dbse, 'h2o2'                  )] = dict(zip(ACTV['%s-%s' % (dbse, 'h2o2')], [+1]))

ACTV['%s-%s'            % (dbse, 'h2o3'                  )] = ['%s-%s-reagent'      % (dbse, 'h2o3')]
RXNM['%s-%s'            % (dbse, 'h2o3'                  )] = dict(zip(ACTV['%s-%s' % (dbse, 'h2o3')], [+1]))

ACTV['%s-%s'            % (dbse, 'h2o4'                  )] = ['%s-%s-reagent'      % (dbse, 'h2o4')]
RXNM['%s-%s'            % (dbse, 'h2o4'                  )] = dict(zip(ACTV['%s-%s' % (dbse, 'h2o4')], [+1]))

ACTV['%s-%s'            % (dbse, 'h2o5'                  )] = ['%s-%s-reagent'      % (dbse, 'h2o5')]
RXNM['%s-%s'            % (dbse, 'h2o5'                  )] = dict(zip(ACTV['%s-%s' % (dbse, 'h2o5')], [+1]))

ACTV['%s-%s'            % (dbse, 'h2o6'                  )] = ['%s-%s-reagent'      % (dbse, 'h2o6')]
RXNM['%s-%s'            % (dbse, 'h2o6'                  )] = dict(zip(ACTV['%s-%s' % (dbse, 'h2o6')], [+1]))

ACTV['%s-%s'            % (dbse, 'Kr'                    )] = ['%s-%s-reagent'      % (dbse, 'Kr')]
RXNM['%s-%s'            % (dbse, 'Kr'                    )] = dict(zip(ACTV['%s-%s' % (dbse, 'Kr')], [+1]))

ACTV['%s-%s'            % (dbse, 'Ne'                    )] = ['%s-%s-reagent'      % (dbse, 'Ne')]
RXNM['%s-%s'            % (dbse, 'Ne'                    )] = dict(zip(ACTV['%s-%s' % (dbse, 'Ne')], [+1]))

ACTV['%s-%s'            % (dbse, 'prism'                 )] = ['%s-%s-reagent'      % (dbse, 'prism')]
RXNM['%s-%s'            % (dbse, 'prism'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'prism')], [+1]))

ACTV['%s-%s'            % (dbse, 'S22_DD_bz_2_pd'        )] = ['%s-%s-reagent'      % (dbse, 'S22_DD_bz_2_pd')]
RXNM['%s-%s'            % (dbse, 'S22_DD_bz_2_pd'        )] = dict(zip(ACTV['%s-%s' % (dbse, 'S22_DD_bz_2_pd')], [+1]))

ACTV['%s-%s'            % (dbse, 'S22_DD_bz_2_t'         )] = ['%s-%s-reagent'      % (dbse, 'S22_DD_bz_2_t')]
RXNM['%s-%s'            % (dbse, 'S22_DD_bz_2_t'         )] = dict(zip(ACTV['%s-%s' % (dbse, 'S22_DD_bz_2_t')], [+1]))

ACTV['%s-%s'            % (dbse, 'S22_DD_bz_me'          )] = ['%s-%s-reagent'      % (dbse, 'S22_DD_bz_me')]
RXNM['%s-%s'            % (dbse, 'S22_DD_bz_me'          )] = dict(zip(ACTV['%s-%s' % (dbse, 'S22_DD_bz_me')], [+1]))

ACTV['%s-%s'            % (dbse, 'S22_DD_ch4_2'          )] = ['%s-%s-reagent'      % (dbse, 'S22_DD_ch4_2')]
RXNM['%s-%s'            % (dbse, 'S22_DD_ch4_2'          )] = dict(zip(ACTV['%s-%s' % (dbse, 'S22_DD_ch4_2')], [+1]))

ACTV['%s-%s'            % (dbse, 'S22_HB_formamide_2'    )] = ['%s-%s-reagent'      % (dbse, 'S22_HB_formamide_2')]
RXNM['%s-%s'            % (dbse, 'S22_HB_formamide_2'    )] = dict(zip(ACTV['%s-%s' % (dbse, 'S22_HB_formamide_2')], [+1]))

ACTV['%s-%s'            % (dbse, 'S22_HB_formic_2'       )] = ['%s-%s-reagent'      % (dbse, 'S22_HB_formic_2')]
RXNM['%s-%s'            % (dbse, 'S22_HB_formic_2'       )] = dict(zip(ACTV['%s-%s' % (dbse, 'S22_HB_formic_2')], [+1]))

ACTV['%s-%s'            % (dbse, 'S22_HB_h2o_2'          )] = ['%s-%s-reagent'      % (dbse, 'S22_HB_h2o_2')]
RXNM['%s-%s'            % (dbse, 'S22_HB_h2o_2'          )] = dict(zip(ACTV['%s-%s' % (dbse, 'S22_HB_h2o_2')], [+1]))

ACTV['%s-%s'            % (dbse, 'S22_HB_nh3_2'          )] = ['%s-%s-reagent'      % (dbse, 'S22_HB_nh3_2')]
RXNM['%s-%s'            % (dbse, 'S22_HB_nh3_2'          )] = dict(zip(ACTV['%s-%s' % (dbse, 'S22_HB_nh3_2')], [+1]))

ACTV['%s-%s'            % (dbse, 'S22_MX_bz_h2o'         )] = ['%s-%s-reagent'      % (dbse, 'S22_MX_bz_h2o')]
RXNM['%s-%s'            % (dbse, 'S22_MX_bz_h2o'         )] = dict(zip(ACTV['%s-%s' % (dbse, 'S22_MX_bz_h2o')], [+1]))

ACTV['%s-%s'            % (dbse, 'S22_MX_bz_hcn'         )] = ['%s-%s-reagent'      % (dbse, 'S22_MX_bz_hcn')]
RXNM['%s-%s'            % (dbse, 'S22_MX_bz_hcn'         )] = dict(zip(ACTV['%s-%s' % (dbse, 'S22_MX_bz_hcn')], [+1]))

ACTV['%s-%s'            % (dbse, 'S22_MX_bz_nh3'         )] = ['%s-%s-reagent'      % (dbse, 'S22_MX_bz_nh3')]
RXNM['%s-%s'            % (dbse, 'S22_MX_bz_nh3'         )] = dict(zip(ACTV['%s-%s' % (dbse, 'S22_MX_bz_nh3')], [+1]))

ACTV['%s-%s'            % (dbse, 'S22_MX_c2h2_c2h4'      )] = ['%s-%s-reagent'      % (dbse, 'S22_MX_c2h2_c2h4')]
RXNM['%s-%s'            % (dbse, 'S22_MX_c2h2_c2h4'      )] = dict(zip(ACTV['%s-%s' % (dbse, 'S22_MX_c2h2_c2h4')], [+1]))

ACTV['%s-%s'            % (dbse, 'thio1'                 )] = ['%s-%s-reagent'      % (dbse, 'thio1')]
RXNM['%s-%s'            % (dbse, 'thio1'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'thio1')], [+1]))

ACTV['%s-%s'            % (dbse, 'thio2'                 )] = ['%s-%s-reagent'      % (dbse, 'thio2')]
RXNM['%s-%s'            % (dbse, 'thio2'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'thio2')], [+1]))

# <<< Reference Values [kcal/mol] >>>
BIND = {}
nan = float('NaN')
BIND['%s-%s'            % (dbse, 'acene1'                )] =      nan
BIND['%s-%s'            % (dbse, 'acene2'                )] =      nan
BIND['%s-%s'            % (dbse, 'alka10'                )] =      nan
BIND['%s-%s'            % (dbse, 'alka2'                 )] =      nan
BIND['%s-%s'            % (dbse, 'alka4'                 )] =      nan
BIND['%s-%s'            % (dbse, 'alka6'                 )] =      nan
BIND['%s-%s'            % (dbse, 'alka8'                 )] =      nan
BIND['%s-%s'            % (dbse, 'alke10'                )] =      nan
BIND['%s-%s'            % (dbse, 'alke2'                 )] =      nan
BIND['%s-%s'            % (dbse, 'alke4'                 )] =      nan
BIND['%s-%s'            % (dbse, 'alke6'                 )] =      nan
BIND['%s-%s'            % (dbse, 'alke8'                 )] =      nan
BIND['%s-%s'            % (dbse, 'Ar'                    )] =      nan
BIND['%s-%s'            % (dbse, 'bag'                   )] =      nan
BIND['%s-%s'            % (dbse, 'boat1'                 )] =      nan
BIND['%s-%s'            % (dbse, 'boat2'                 )] =      nan
BIND['%s-%s'            % (dbse, 'book1'                 )] =      nan
BIND['%s-%s'            % (dbse, 'book2'                 )] =      nan
BIND['%s-%s'            % (dbse, 'cage'                  )] =      nan
BIND['%s-%s'            % (dbse, 'cyclic'                )] =      nan
BIND['%s-%s'            % (dbse, 'h2o'                   )] =      nan
BIND['%s-%s'            % (dbse, 'h2o2'                  )] =      nan
BIND['%s-%s'            % (dbse, 'h2o3'                  )] =      nan
BIND['%s-%s'            % (dbse, 'h2o4'                  )] =      nan
BIND['%s-%s'            % (dbse, 'h2o5'                  )] =      nan
BIND['%s-%s'            % (dbse, 'h2o6'                  )] =      nan
BIND['%s-%s'            % (dbse, 'Kr'                    )] =      nan
BIND['%s-%s'            % (dbse, 'Ne'                    )] =      nan
BIND['%s-%s'            % (dbse, 'prism'                 )] =      nan
BIND['%s-%s'            % (dbse, 'S22_DD_bz_2_pd'        )] =      nan
BIND['%s-%s'            % (dbse, 'S22_DD_bz_2_t'         )] =      nan
BIND['%s-%s'            % (dbse, 'S22_DD_bz_me'          )] =      nan
BIND['%s-%s'            % (dbse, 'S22_DD_ch4_2'          )] =      nan
BIND['%s-%s'            % (dbse, 'S22_HB_formamide_2'    )] =      nan
BIND['%s-%s'            % (dbse, 'S22_HB_formic_2'       )] =      nan
BIND['%s-%s'            % (dbse, 'S22_HB_h2o_2'          )] =      nan
BIND['%s-%s'            % (dbse, 'S22_HB_nh3_2'          )] =      nan
BIND['%s-%s'            % (dbse, 'S22_MX_bz_h2o'         )] =      nan
BIND['%s-%s'            % (dbse, 'S22_MX_bz_hcn'         )] =      nan
BIND['%s-%s'            % (dbse, 'S22_MX_bz_nh3'         )] =      nan
BIND['%s-%s'            % (dbse, 'S22_MX_c2h2_c2h4'      )] =      nan
BIND['%s-%s'            % (dbse, 'thio1'                 )] =      nan
BIND['%s-%s'            % (dbse, 'thio2'                 )] =      nan

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, 'acene1'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'acene1'                )] = """ """
TAGL['%s-%s'            % (dbse, 'acene2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'acene2'                )] = """ """
TAGL['%s-%s'            % (dbse, 'alka10'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alka10'                )] = """ """
TAGL['%s-%s'            % (dbse, 'alka2'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alka2'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'alka4'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alka4'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'alka6'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alka6'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'alka8'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alka8'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'alke10'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alke10'                )] = """ """
TAGL['%s-%s'            % (dbse, 'alke2'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alke2'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'alke4'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alke4'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'alke6'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alke6'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'alke8'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alke8'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'Ar'                    )] = """-1 0 """
TAGL['%s-%s-reagent'    % (dbse, 'Ar'                    )] = """-1 0 """
TAGL['%s-%s'            % (dbse, 'bag'                   )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'bag'                   )] = """ """
TAGL['%s-%s'            % (dbse, 'boat1'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'boat1'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'boat2'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'boat2'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'book1'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'book1'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'book2'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'book2'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'cage'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'cage'                  )] = """ """
TAGL['%s-%s'            % (dbse, 'cyclic'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'cyclic'                )] = """ """
TAGL['%s-%s'            % (dbse, 'h2o'                   )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'h2o'                   )] = """ """
TAGL['%s-%s'            % (dbse, 'h2o2'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'h2o2'                  )] = """ """
TAGL['%s-%s'            % (dbse, 'h2o3'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'h2o3'                  )] = """ """
TAGL['%s-%s'            % (dbse, 'h2o4'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'h2o4'                  )] = """ """
TAGL['%s-%s'            % (dbse, 'h2o5'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'h2o5'                  )] = """ """
TAGL['%s-%s'            % (dbse, 'h2o6'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'h2o6'                  )] = """ """
TAGL['%s-%s'            % (dbse, 'Kr'                    )] = """-1 0 """
TAGL['%s-%s-reagent'    % (dbse, 'Kr'                    )] = """-1 0 """
TAGL['%s-%s'            % (dbse, 'Ne'                    )] = """-1 0 """
TAGL['%s-%s-reagent'    % (dbse, 'Ne'                    )] = """-1 0 """
TAGL['%s-%s'            % (dbse, 'prism'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'prism'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'S22_DD_bz_2_pd'        )] = """0 1 set S22 index 11 DD-4 benzene_dimer c2h """
TAGL['%s-%s-reagent'    % (dbse, 'S22_DD_bz_2_pd'        )] = """0 1 set S22 index 11 DD-4 benzene_dimer c2h """
TAGL['%s-%s'            % (dbse, 'S22_DD_bz_2_t'         )] = """0 1 set S22 index 20 MX-5 benzene_dimer c2v """
TAGL['%s-%s-reagent'    % (dbse, 'S22_DD_bz_2_t'         )] = """0 1 set S22 index 20 MX-5 benzene_dimer c2v """
TAGL['%s-%s'            % (dbse, 'S22_DD_bz_me'          )] = """0 1 set S22 index 10 DD-3 benzene_methane c3 """
TAGL['%s-%s-reagent'    % (dbse, 'S22_DD_bz_me'          )] = """0 1 set S22 index 10 DD-3 benzene_methane c3 """
TAGL['%s-%s'            % (dbse, 'S22_DD_ch4_2'          )] = """0 1 set S22 index 8 DD-1 methane_dimer d3d """
TAGL['%s-%s-reagent'    % (dbse, 'S22_DD_ch4_2'          )] = """0 1 set S22 index 8 DD-1 methane_dimer d3d """
TAGL['%s-%s'            % (dbse, 'S22_HB_formamide_2'    )] = """0 1 set S22 index 4 HB-4 formamide_dimer c2h """
TAGL['%s-%s-reagent'    % (dbse, 'S22_HB_formamide_2'    )] = """0 1 set S22 index 4 HB-4 formamide_dimer c2h """
TAGL['%s-%s'            % (dbse, 'S22_HB_formic_2'       )] = """0 1 set S22 index 3 HB-3 formic_acid_dimer c2h """
TAGL['%s-%s-reagent'    % (dbse, 'S22_HB_formic_2'       )] = """0 1 set S22 index 3 HB-3 formic_acid_dimer c2h """
TAGL['%s-%s'            % (dbse, 'S22_HB_h2o_2'          )] = """0 1 set S22 index 2 HB-2 water_dimer cs """
TAGL['%s-%s-reagent'    % (dbse, 'S22_HB_h2o_2'          )] = """0 1 set S22 index 2 HB-2 water_dimer cs """
TAGL['%s-%s'            % (dbse, 'S22_HB_nh3_2'          )] = """0 1 set S22 index 1 HB-1 ammonia_dimer c2h """
TAGL['%s-%s-reagent'    % (dbse, 'S22_HB_nh3_2'          )] = """0 1 set S22 index 1 HB-1 ammonia_dimer c2h """
TAGL['%s-%s'            % (dbse, 'S22_MX_bz_h2o'         )] = """0 1 set S22 index 17 MX-2 benzene_water cs """
TAGL['%s-%s-reagent'    % (dbse, 'S22_MX_bz_h2o'         )] = """0 1 set S22 index 17 MX-2 benzene_water cs """
TAGL['%s-%s'            % (dbse, 'S22_MX_bz_hcn'         )] = """0 1 set S22 index 19 MX-4 benzene_hcn cs """
TAGL['%s-%s-reagent'    % (dbse, 'S22_MX_bz_hcn'         )] = """0 1 set S22 index 19 MX-4 benzene_hcn cs """
TAGL['%s-%s'            % (dbse, 'S22_MX_bz_nh3'         )] = """0 1 set S22 index 18 MX-3 benzene_ammonia cs """
TAGL['%s-%s-reagent'    % (dbse, 'S22_MX_bz_nh3'         )] = """0 1 set S22 index 18 MX-3 benzene_ammonia cs """
TAGL['%s-%s'            % (dbse, 'S22_MX_c2h2_c2h4'      )] = """0 1 set S22 index 16 MX-1 ethene_ethine c2v """
TAGL['%s-%s-reagent'    % (dbse, 'S22_MX_c2h2_c2h4'      )] = """0 1 set S22 index 16 MX-1 ethene_ethine c2v """
TAGL['%s-%s'            % (dbse, 'thio1'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'thio1'                 )] = """ """
TAGL['%s-%s'            % (dbse, 'thio2'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'thio2'                 )] = """ """

# <<< Molecule Specifications >>>
BENCH12_acene1 = input.process_input("""
molecule dimer {
0 1
C       -2.21210099    -1.64058681     0.00000000
C       -0.81694099    -1.64058681     0.00000000
C       -0.11940299    -0.43283581     0.00000000
C       -0.81705699     0.77567319    -0.00119900
C       -2.21188199     0.77559519    -0.00167800
C       -2.90948299    -0.43261081    -0.00068200
H       -2.76185999    -2.59290381     0.00045000
H       -0.26743299    -2.59309981     0.00131500
H        0.98027701    -0.43275581     0.00063400
H       -0.26685699     1.72781619    -0.00125800
H       -2.76200399     1.72787619    -0.00263100
H       -4.00908699    -0.43242781    -0.00086200
units angstrom
}
""", 0)

BENCH12_acene2 = input.process_input("""
molecule dimer {
0 1
C        2.42369900     0.70669500     0.00000000
C       -2.42369900    -0.70669500     0.00000000
C        2.42369900    -0.70669500     0.00000000
C       -2.42369900     0.70669500     0.00000000
C        1.24011400     1.39713900     0.00000000
C        1.24011400    -1.39713900     0.00000000
C       -1.24011400     1.39713900     0.00000000
C       -1.24011400    -1.39713900     0.00000000
C        0.00000000     0.71391600     0.00000000
C        0.00000000    -0.71391600     0.00000000
H       -1.24044800     2.47918000     0.00000000
H        3.36402200     1.24175400     0.00000000
H        1.23857800     2.48008100     0.00000000
H        1.24038500    -2.47913700     0.00000000
H        3.36491600    -1.23926600     0.00000000
H       -1.23863300    -2.48004900     0.00000000
H       -3.36410000    -1.24169900     0.00000000
H       -3.36498700     1.23929300     0.00000000
units angstrom
}
""", 0)

BENCH12_alka10 = input.process_input("""
molecule dimer {
0 1
C        0.40267235     0.64977904     0.00000000
H        1.06210212     0.66547832     0.87441225
H        1.06210212     0.66547832    -0.87441225
C       -0.46431471     1.90898497     0.00000000
H       -1.12377520     1.89323444    -0.87442948
H       -1.12377520     1.89323444     0.87442948
C        0.34027906     3.20874936     0.00000000
H        0.99988197     3.22611504    -0.87444682
H        0.99988197     3.22611504     0.87444682
C       -0.52663804     4.46801317     0.00000000
H       -1.18529556     4.45181950     0.87386445
H       -1.18529556     4.45181950    -0.87386445
C        0.28501763     5.76202328     0.00000000
H       -0.36242574     6.64019737     0.00000000
H        0.92798572     5.82534568     0.88054822
H        0.92798572     5.82534568    -0.88054822
C       -0.40267235    -0.64977904     0.00000000
H       -1.06210212    -0.66547832     0.87441225
H       -1.06210212    -0.66547832    -0.87441225
C        0.46431471    -1.90898497     0.00000000
H        1.12377520    -1.89323444     0.87442948
H        1.12377520    -1.89323444    -0.87442948
C       -0.34027906    -3.20874936     0.00000000
H       -0.99988197    -3.22611504    -0.87444682
H       -0.99988197    -3.22611504     0.87444682
C        0.52663804    -4.46801317     0.00000000
H        1.18529556    -4.45181950    -0.87386445
H        1.18529556    -4.45181950     0.87386445
C       -0.28501763    -5.76202328     0.00000000
H        0.36242574    -6.64019737     0.00000000
H       -0.92798572    -5.82534568    -0.88054822
H       -0.92798572    -5.82534568     0.88054822
units angstrom
}
""", 0)

BENCH12_alka2 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -0.76370350
H        0.00000000     1.01617586    -1.16095589
H        0.88003411    -0.50808793    -1.16095589
H       -0.88003411    -0.50808793    -1.16095589
C        0.00000000     0.00000000     0.76370350
H        0.88003411     0.50808793     1.16095589
H       -0.88003411     0.50808793     1.16095589
H       -0.00000000    -1.01617586     1.16095589
units angstrom
}
""", 0)

BENCH12_alka4 = input.process_input("""
molecule dimer {
0 1
C        0.51479500    -0.56820117     0.00000000
H        1.17232619    -0.45994153     0.88262758
H        1.17232619    -0.45994153    -0.88262758
C       -0.12005863    -1.96051588     0.00000000
H        0.64293252    -2.75572018     0.00000000
H       -0.75694888    -2.10907685     0.88899667
H       -0.75694888    -2.10907685    -0.88899667
C       -0.51479500     0.56820117     0.00000000
H       -1.17232619     0.45994153    -0.88262758
H       -1.17232619     0.45994153     0.88262758
C        0.12005863     1.96051588     0.00000000
H       -0.64293252     2.75572018     0.00000000
H        0.75694888     2.10907685    -0.88899667
H        0.75694888     2.10907685     0.88899667
units angstrom
}
""", 0)

BENCH12_alka6 = input.process_input("""
molecule dimer {
0 1
C        0.21092141     3.21935408     0.00000000
C       -0.55072178     1.89187840     0.00000000
H       -1.21588133     1.84802196     0.88164500
H       -0.47592852     4.08118449     0.00000000
H        0.85892441     3.30920214     0.88746775
H        0.85892441     3.30920214    -0.88746775
C        0.36718989     0.66476129     0.00000000
H       -1.21588133     1.84802196    -0.88164500
H        1.03339427     0.70846132    -0.88345369
H        1.03339427     0.70846132     0.88345369
C       -0.36718989    -0.66476129     0.00000000
C        0.55072178    -1.89187840     0.00000000
H       -1.03339427    -0.70846132    -0.88345369
H       -1.03339427    -0.70846132     0.88345369
C       -0.21092141    -3.21935408     0.00000000
H        1.21588133    -1.84802196     0.88164500
H        1.21588133    -1.84802196    -0.88164500
H       -0.85892441    -3.30920214     0.88746775
H       -0.85892441    -3.30920214    -0.88746775
H        0.47592852    -4.08118449     0.00000000
units angstrom
}
""", 0)

BENCH12_alka8 = input.process_input("""
molecule dimer {
0 1
C        0.25853799     4.49476029     0.00000000
C       -0.54104699     3.18978676     0.00000000
H       -1.20719420     3.16511861     0.88164500
H       -0.40318847     5.37602789     0.00000000
H        0.90886126     4.56589543     0.88746775
H        0.90886126     4.56589543    -0.88746775
C        0.34111758     1.93672496     0.00000000
H       -1.20719420     3.16511861    -0.88164500
H        1.00830467     1.96120666    -0.88345369
H        1.00830467     1.96120666     0.88345369
C       -0.44965914     0.62446618     0.00000000
C        0.44965914    -0.62446618     0.00000000
H       -1.11728738     0.60043427    -0.88238195
H       -1.11728738     0.60043427     0.88238195
C       -0.34111758    -1.93672496     0.00000000
H        1.11728738    -0.60043427     0.88238195
H        1.11728738    -0.60043427    -0.88238195
H       -1.00830467    -1.96120666     0.88345369
H       -1.00830467    -1.96120666    -0.88345369
C        0.54104699    -3.18978676     0.00000000
C       -0.25853799    -4.49476029     0.00000000
H        1.20719420    -3.16511861     0.88164500
H        1.20719420    -3.16511861    -0.88164500
H       -0.90886126    -4.56589543    -0.88746775
H       -0.90886126    -4.56589543     0.88746775
H        0.40318847    -5.37602789     0.00000000
units angstrom
}
""", 0)

BENCH12_alke10 = input.process_input("""
molecule dimer {
0 1
H       -0.03635200     0.00000000    -0.09668700
C       -0.05408000     0.00000000     0.99258700
C        1.08386800     0.00000000     1.71595900
H       -1.03434700    -0.00000100     1.47263100
H        2.04457700     0.00000000     1.19137600
C        1.13596400     0.00000000     3.16322500
C        2.28721200     0.00000000     3.88723800
H        0.17781700     0.00000000     3.69317200
H        3.24225100     0.00000000     3.35159300
C        2.34560400     0.00000000     5.32584600
C        3.50006700     0.00000000     6.04915400
H        1.39086200     0.00000000     5.86174000
H        4.45480900     0.00000000     5.51326000
C        3.55845900     0.00000000     7.48776200
C        4.70970800     0.00000000     8.21177500
H        2.60342100     0.00000000     8.02340700
H        5.66785500     0.00000000     7.68182800
C        4.76180400     0.00000000     9.65904100
C        5.89975200     0.00000000    10.38241300
H        3.80109500     0.00000000    10.18362400
H        6.88001800     0.00000000     9.90236900
H        5.88202300    -0.00000100    11.47168700
units angstrom
}
""", 0)

BENCH12_alke2 = input.process_input("""
molecule dimer {
0 1
C        0.00000000    -0.66757800    -2.12465900
C        0.00000000     0.66757800    -2.12465900
H        0.92362100    -1.23225300    -2.12618500
H       -0.92362100    -1.23225300    -2.12618500
H       -0.92362100     1.23225300    -2.12618500
H        0.92362100     1.23225300    -2.12618500
units angstrom
}
""", 0)

BENCH12_alke4 = input.process_input("""
molecule dimer {
0 1
H       -0.02710600     0.00000000    -0.01429600
C       -0.01998300     0.00000000     1.07543800
C        1.12986600     0.00000000     1.77174400
H       -0.98974000     0.00000000     1.57688300
H        2.08205900     0.00000000     1.23175800
C        1.20840200     0.00000000     3.22825600
C        2.35825100     0.00000000     3.92456200
H        0.25621000     0.00000000     3.76824200
H        3.32800800     0.00000000     3.42311700
H        2.36537500     0.00000000     5.01429600
units angstrom
}
""", 0)

BENCH12_alke6 = input.process_input("""
molecule dimer {
0 1
H       -0.02920300     0.00000000    -0.04214200
C       -0.03456600     0.00000000     1.04735400
C        1.11009000     0.00000000     1.75717000
H       -1.00959700     0.00000000     1.53813600
H        2.06568600     0.00000000     1.22342900
C        1.17692500     0.00000000     3.20692100
C        2.33047800     0.00000000     3.91807900
H        0.22335700     0.00000000     3.74536600
H        3.28404600     0.00000000     3.37963400
C        2.39731300     0.00000000     5.36783000
C        3.54196900     0.00000000     6.07764600
H        1.44171700     0.00000000     5.90157200
H        4.51700000     0.00000000     5.58686400
H        3.53660600     0.00000000     7.16714200
units angstrom
}
""", 0)

BENCH12_alke8 = input.process_input("""
molecule dimer {
0 1
H       -0.03267500     0.00000000    -0.06940800
C       -0.04516000     0.00000000     1.01994800
C        1.09607000     0.00000000     1.73733100
H       -1.02260300     0.00000000     1.50575100
H        2.05459000     0.00000000     1.20870600
C        1.15455000     0.00000000     3.18528400
C        2.30705400     0.00000000     3.90435100
H        0.19816000     0.00000000     3.71860300
H        3.26106800     0.00000000     3.36691600
C        2.36948400     0.00000000     5.34564900
C        3.52198700     0.00000000     6.06471600
H        1.41546900     0.00000000     5.88308300
H        4.47837800     0.00000000     5.53139700
C        3.58046700     0.00000000     7.51266900
C        4.72169700     0.00000000     8.23005200
H        2.62194700     0.00000000     8.04129400
H        5.69914000     0.00000000     7.74424900
H        4.70921200     0.00000000     9.31940800
units angstrom
}
""", 0)

BENCH12_Ar = input.process_input("""
molecule dimer {
0 1
Ar       0.00000000     0.00000000     0.00000000
units angstrom
}
""", 0)

BENCH12_bag = input.process_input("""
molecule dimer {
0 1
O       -1.02199800     1.60831800    -0.45667800
H       -0.03478700     1.69574200    -0.37604200
H       -1.32678900     2.46174400    -0.77816300
O       -0.09409400    -1.92801900     0.49396100
H       -0.53052100    -1.35219800     1.14518000
H       -0.58296100    -1.74114400    -0.32694300
O        1.62181300     1.69683700    -0.20757000
H        1.95015200     0.76515400    -0.16297800
H        1.96618900     2.11133200     0.58869800
O        2.34085300    -0.86771800     0.06047900
H        1.48693800    -1.33684800     0.23150500
H        2.78403500    -1.37920100    -0.62125300
O       -1.33564800     0.18185300     1.95656400
H       -1.40191400     0.79271400     1.20173800
H       -2.19955700     0.19488200     2.37668200
O       -1.33146000    -0.85096500    -1.80507500
H       -1.32899000     0.06249800    -1.46551700
H       -2.21752000    -0.99712800    -2.14635300
units angstrom
}
""", 0)

BENCH12_boat1 = input.process_input("""
molecule dimer {
0 1
O       -1.72674400    -1.94937300    -0.28696700
H       -2.04961100    -1.02159300    -0.32910600
H       -2.49471900    -2.47594700    -0.05241700
O       -2.59026600     0.62259900    -0.36163600
H       -1.94834000     1.27648400    -0.00529500
H       -2.87658700     0.99098200    -1.20154400
O       -0.79461000     2.41504100     0.61061700
H        0.13441000     2.26696700     0.32798500
H       -0.75442700     2.53636700     1.56255200
O        2.59016300    -0.62273600    -0.36190100
H        2.87625100    -0.99116700    -1.20186600
H        1.94826400    -1.27655500    -0.00538100
O        0.79455700    -2.41484700     0.61098900
H       -0.13445100    -2.26674200     0.32830100
H        0.75437800    -2.53581200     1.56297100
O        1.72692200     1.94930100    -0.28670600
H        2.04959200     1.02146900    -0.32914400
H        2.49505400     2.47567300    -0.05222900
units angstrom
}
""", 0)

BENCH12_boat2 = input.process_input("""
molecule dimer {
0 1
O       -1.38674411    -2.13681196    -0.37562259
H       -0.48634999    -2.17895521     0.01536018
H       -1.77328716    -3.00242730    -0.22364016
O       -1.13219127     2.21240776     0.64824427
H       -1.29929148     2.23863212     1.59379084
H       -1.72071662     1.50204924     0.30641457
O       -2.72623407     0.22572348    -0.29552824
H       -2.27361904    -0.64532669    -0.32433470
H       -3.06509937     0.35947845    -1.18468601
O        1.38674221     2.13681088    -0.37562863
H        1.77329342     3.00242205    -0.22364361
H        0.48633699     2.17897063     0.01532572
O        1.13215726    -2.21237218     0.64832349
H        1.29922122    -2.23853109     1.59387834
H        1.72069165    -1.50203363     0.30646786
O        2.72621384    -0.22573640    -0.29553397
H        3.06503646    -0.35951676    -1.18470419
H        2.27360404     0.64531609    -0.32434222
units angstrom
}
""", 0)

BENCH12_book1 = input.process_input("""
molecule dimer {
0 1
O        0.11761300     1.55200200     0.86493100
H        0.96325100     1.52043000     0.35419900
H        0.27939600     2.15144900     1.59904500
O       -0.02656300    -1.36928700     0.91920900
H       -0.04796800    -0.42097200     1.11892700
H       -0.87592000    -1.52924100     0.46385900
O       -2.33111700     1.36041800    -0.43366500
H       -1.46601800     1.54388000    -0.01880800
H       -2.34390900     1.88470800    -1.23876600
O       -2.47246500    -1.40646100    -0.43402900
H       -2.54678700    -0.43378600    -0.50026700
H       -3.28021400    -1.69060900     0.00214200
O        2.41910100     1.20412500    -0.46668500
H        2.49230900     0.22056700    -0.49973000
H        2.55739900     1.49934200    -1.37062000
O        2.28641700    -1.48073200    -0.44513000
H        1.42583000    -1.57095100     0.03130600
H        2.89873600    -2.05534300     0.02166200
units angstrom
}
""", 0)

BENCH12_book2 = input.process_input("""
molecule dimer {
0 1
O        2.40320534     1.32706594    -0.42312413
H        3.28247053     1.52755300    -0.09144602
H        2.40101026     0.36107264    -0.57991790
O        2.06820382    -1.40350163    -0.65215998
H        1.86446942    -1.80746273    -1.50019296
H        1.28096132    -1.56788983    -0.09729097
O       -0.16435795    -1.55958995     1.00265246
H       -1.01292419    -1.52303353     0.49858555
H       -0.28266618    -2.23257352     1.67874174
O       -2.31816171    -1.15518688    -0.54535437
H       -3.22909379    -1.38580965    -0.34437983
H       -2.30926544    -0.17053375    -0.61709369
O       -2.06411969     1.52008440    -0.49855150
H       -1.85291346     2.04318067    -1.27644018
H       -1.26825985     1.58272526     0.08372064
O        0.08217709     1.35119088     1.13167187
H        0.05811694     0.40256100     1.32906529
H        0.92799546     1.46558737     0.65660106
units angstrom
}
""", 0)

BENCH12_cage = input.process_input("""
molecule dimer {
0 1
O        0.86429700    -1.71381800    -0.47336700
H        1.68527300    -1.20159100    -0.30198200
H        1.14997000    -2.58689700    -0.75629900
O       -0.73069700     0.40566900    -1.63668400
H       -0.26478000    -0.43576700    -1.50730100
H       -1.61069500     0.24680500    -1.24492900
O        0.59012400     1.73858300     0.25043100
H        0.06749400     1.37261900    -0.51268400
H        0.49869400     2.69429100     0.20132100
O       -0.70229900    -0.49472100     1.65945100
H       -0.20123800    -1.09653500     1.08595400
H       -0.28260900     0.36322900     1.49278700
O        2.79584000     0.11344200     0.07992500
H        2.14917800     0.84175600     0.18112500
H        3.29643200     0.10246400     0.90040600
O       -2.89018000     0.01527900     0.06304800
H       -2.25688400    -0.24808500     0.76309500
H       -3.64752400    -0.56775600     0.15607900
units angstrom
}
""", 0)

BENCH12_cyclic = input.process_input("""
molecule dimer {
0 1
O        0.00000000     2.69654700     0.13950200
H       -0.12861000     3.23017500     0.92778900
H       -0.84251400     2.20238000     0.02595300
O       -2.33527800    -1.34827400     0.13950200
H       -1.48606000    -1.83082800     0.02595300
H       -2.73310800    -1.72646700     0.92778900
O        2.33527800     1.34827400    -0.13950200
H        2.73310800     1.72646700    -0.92778900
H        1.48606000     1.83082800    -0.02595300
O        0.00000000    -2.69654700    -0.13950200
H        0.84251400    -2.20238000    -0.02595300
H        0.12861000    -3.23017500    -0.92778900
O        2.33527800    -1.34827400     0.13950200
H        2.86171800    -1.50370800     0.92778900
H        2.32857400    -0.37155100     0.02595300
O       -2.33527800     1.34827400    -0.13950200
H       -2.86171800     1.50370800    -0.92778900
H       -2.32857400     0.37155100    -0.02595300
units angstrom
}
""", 0)

BENCH12_h2o = input.process_input("""
molecule dimer {
0 1
O        0.00000000     0.00000000    -0.12789657
H        0.00000000    -1.42990030     1.01490567
H        0.00000000     1.42990030     1.01490567
units angstrom
}
""", 0)

BENCH12_h2o2 = input.process_input("""
molecule dimer {
0 1
O       -1.55100700    -0.11452000     0.00000000
H       -1.93425900     0.76250300     0.00000000
H       -0.59967700     0.04071200     0.00000000
O        1.35062500     0.11146900     0.00000000
H        1.68039800    -0.37374100    -0.75856100
H        1.68039800    -0.37374100     0.75856100
units angstrom
}
""", 0)

BENCH12_h2o3 = input.process_input("""
molecule dimer {
0 1
O       -1.29527542    -0.95244604    -0.08439711
H       -1.96684629    -1.20450472     0.55440905
H       -1.21891600     0.01667075    -0.01273795
O       -0.17970784     1.59357918     0.12071627
H       -0.09171077     2.24303363    -0.58191728
H        0.62465016     1.04596562     0.06646441
O        1.47586351    -0.63755538    -0.07513260
H        2.03414372    -1.07560919     0.57227133
H        0.61081072    -1.08144807    -0.01749750
units angstrom
}
""", 0)

BENCH12_h2o4 = input.process_input("""
molecule dimer {
0 1
O        1.94056518    -0.00328371     0.03503401
H        1.37319186     0.79405277    -0.01037782
H        2.67990412     0.16977815    -0.55264726
O        0.00328371     1.94056518     0.03503401
H       -0.79405277     1.37319186    -0.01037782
H       -0.16977815     2.67990412    -0.55264726
O       -1.94056518     0.00328371     0.03503401
H       -1.37319186    -0.79405277    -0.01037782
H       -2.67990412    -0.16977815    -0.55264726
O       -0.00328371    -1.94056518     0.03503401
H        0.79405277    -1.37319186    -0.01037782
H        0.16977815    -2.67990412    -0.55264726
units angstrom
}
""", 0)

BENCH12_h2o5 = input.process_input("""
molecule dimer {
0 1
O        0.71936714     2.21304547    -0.12177106
H       -0.23169873     1.97440618    -0.09174602
H        0.80185210     3.00522157     0.41416587
O        2.32702474    -0.00028534    -0.12177106
H        1.80616977     0.83048845    -0.09174602
H        3.10591824     0.16606268     0.41416587
O        0.71881557    -2.21321542    -0.12177106
H        1.34797537    -1.46112968    -0.09174602
H        1.11771327    -2.90258279     0.41416587
O       -1.88276996    -1.36755061    -0.12177106
H       -0.97307285    -1.73350986    -0.09174602
H       -2.41513112    -1.95995109     0.41416587
O       -1.88242907     1.36802906    -0.12177106
H       -1.94936513     0.38976808    -0.09174602
H       -2.61034406     1.69127280     0.41416587
units angstrom
}
""", 0)

BENCH12_h2o6 = input.process_input("""
molecule dimer {
0 1
O        0.00000000     2.69654700     0.13950200
H       -0.12861000     3.23017500     0.92778900
H       -0.84251400     2.20238000     0.02595300
O       -2.33527800    -1.34827400     0.13950200
H       -1.48606000    -1.83082800     0.02595300
H       -2.73310800    -1.72646700     0.92778900
O        2.33527800     1.34827400    -0.13950200
H        2.73310800     1.72646700    -0.92778900
H        1.48606000     1.83082800    -0.02595300
O        0.00000000    -2.69654700    -0.13950200
H        0.84251400    -2.20238000    -0.02595300
H        0.12861000    -3.23017500    -0.92778900
O        2.33527800    -1.34827400     0.13950200
H        2.86171800    -1.50370800     0.92778900
H        2.32857400    -0.37155100     0.02595300
O       -2.33527800     1.34827400    -0.13950200
H       -2.86171800     1.50370800    -0.92778900
H       -2.32857400     0.37155100    -0.02595300
units angstrom
}
""", 0)

BENCH12_Kr = input.process_input("""
molecule dimer {
0 1
Kr       0.00000000     0.00000000     0.00000000
units angstrom
}
""", 0)

BENCH12_Ne = input.process_input("""
molecule dimer {
0 1
Ne       0.00000000     0.00000000     0.00000000
units angstrom
}
""", 0)

BENCH12_prism = input.process_input("""
molecule dimer {
0 1
O       -1.50216900    -0.19135900     1.43492700
H       -0.60105400    -0.59697200     1.55371800
H       -2.00669800    -0.42232700     2.21984700
O       -1.74457500    -0.38234800    -1.30914400
H       -1.88894100    -0.47965300    -0.34762400
H       -2.51683500    -0.76676500    -1.73376600
O       -0.56040900     2.01783000    -0.12198400
H       -0.94772000     1.53356700     0.62522800
H       -0.98983100     1.59273600    -0.87741900
O        0.96480300    -1.16576500     1.43998700
H        0.97955700    -1.52204100     0.52783300
H        1.54222400    -0.39369200     1.34437300
O        0.97470500    -1.40150300    -1.33597000
H        0.06516100    -1.11895100    -1.52288600
H        1.47070900    -0.57093300    -1.27771000
O        2.00228000     1.05782400    -0.12450200
H        1.14163700     1.53226600    -0.14012100
H        2.67471600     1.73534200    -0.23799500
units angstrom
}
""", 0)

BENCH12_S22_DD_bz_2_pd = input.process_input("""
molecule dimer {
0 1
C       -1.04782520    -1.42167360     0.00000000
C       -1.45450340    -0.85544590     1.20620480
C       -1.45450340    -0.85544590    -1.20620480
C       -2.26679700     0.27716100     1.20695390
C       -2.67147810     0.84502110     0.00000000
C       -2.26679700     0.27716100    -1.20695390
H       -1.13385340    -1.29205930    -2.14231500
H       -2.58249430     0.71630660    -2.14379770
H       -3.30304220     1.72327000     0.00000000
H       -2.58249430     0.71630660     2.14379770
H       -1.13385340    -1.29205930     2.14231500
H       -0.40602530    -2.29190490     0.00000000
C        1.04782520     1.42167360     0.00000000
C        1.45450340     0.85544590    -1.20620480
C        1.45450340     0.85544590     1.20620480
C        2.26679700    -0.27716100    -1.20695390
C        2.67147810    -0.84502110     0.00000000
C        2.26679700    -0.27716100     1.20695390
H        0.40602530     2.29190490     0.00000000
H        1.13385340     1.29205930     2.14231500
H        2.58249430    -0.71630660     2.14379770
H        3.30304220    -1.72327000     0.00000000
H        2.58249430    -0.71630660    -2.14379770
H        1.13385340     1.29205930    -2.14231500
units angstrom
}
""", 0)

BENCH12_S22_DD_bz_2_t = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     1.05903530
C        0.00000000    -1.20600840     1.75767420
C        0.00000000    -1.20717670     3.15159050
C        0.00000000     0.00000000     3.84857510
C        0.00000000     1.20717670     3.15159050
C        0.00000000     1.20600840     1.75767420
H        0.00000000     0.00000000    -0.02158050
H        0.00000000    -2.14163870     1.21442170
H        0.00000000    -2.14356570     3.69299530
H        0.00000000     0.00000000     4.93014990
H        0.00000000     2.14356570     3.69299530
H        0.00000000     2.14163870     1.21442170
C       -1.39406330     0.00000000    -2.45415240
C       -0.69704680     1.20723780    -2.45462770
C        0.69704680     1.20723780    -2.45462770
C        1.39406330     0.00000000    -2.45415240
C        0.69704680    -1.20723780    -2.45462770
C       -0.69704680    -1.20723780    -2.45462770
H       -2.47539950     0.00000000    -2.45032210
H       -1.23823210     2.14356550    -2.45367640
H        1.23823210     2.14356550    -2.45367640
H        2.47539950     0.00000000    -2.45032210
H        1.23823210    -2.14356550    -2.45367640
H       -1.23823210    -2.14356550    -2.45367640
units angstrom
}
""", 0)

BENCH12_S22_DD_bz_me = input.process_input("""
molecule dimer {
0 1
C        1.39321780     0.03629130    -0.63328030
C        0.72803640    -1.18840150    -0.63330170
C       -0.66517970    -1.22470770    -0.63328030
C       -1.39320410    -0.03629720    -0.63330170
C       -0.72803810     1.18841630    -0.63328030
C        0.66516770     1.22469870    -0.63330170
H        2.47427370     0.06444840    -0.63172400
H        1.29295880    -2.11054090    -0.63174010
H       -1.18132290    -2.17500810    -0.63172400
H       -2.47426140    -0.06446470    -0.63174010
H       -1.29295080     2.11055960    -0.63172400
H        1.18130260     2.17500560    -0.63174010
C        0.00000000     0.00000000     3.08261950
H        0.58687760     0.83817420     3.44637720
H       -1.01931890     0.08916380     3.44637720
H        0.00000000     0.00000000     1.99666970
H        0.43244130    -0.92733800     3.44637720
units angstrom
}
""", 0)

BENCH12_S22_DD_ch4_2 = input.process_input("""
molecule dimer {
0 1
C        0.00000000    -0.00014000     1.85916100
H       -0.88855100     0.51306000     1.49468500
H        0.88855100     0.51306000     1.49468500
H        0.00000000    -1.02633900     1.49486800
H        0.00000000     0.00008900     2.94828400
C        0.00000000     0.00014000    -1.85916100
H        0.00000000    -0.00008900    -2.94828400
H       -0.88855100    -0.51306000    -1.49468500
H        0.88855100    -0.51306000    -1.49468500
H        0.00000000     1.02633900    -1.49486800
units angstrom
}
""", 0)

BENCH12_S22_HB_formamide_2 = input.process_input("""
molecule dimer {
0 1
C       -2.01864900     0.05288300     0.00000000
O       -1.45220000     1.14363400     0.00000000
N       -1.40777000    -1.14248400     0.00000000
H       -1.96459600    -1.97703600     0.00000000
H       -0.38724400    -1.20778200     0.00000000
H       -3.11706100    -0.01370100     0.00000000
C        2.01864900    -0.05288300     0.00000000
O        1.45220000    -1.14363400     0.00000000
N        1.40777000     1.14248400     0.00000000
H        1.96459600     1.97703600     0.00000000
H        0.38724400     1.20778200     0.00000000
H        3.11706100     0.01370100     0.00000000
units angstrom
}
""", 0)

BENCH12_S22_HB_formic_2 = input.process_input("""
molecule dimer {
0 1
C       -1.88889600    -0.17969200     0.00000000
O       -1.49328000     1.07368900     0.00000000
O       -1.17043500    -1.16659000     0.00000000
H       -2.97948800    -0.25882900     0.00000000
H       -0.49883300     1.10719500     0.00000000
C        1.88889600     0.17969200     0.00000000
O        1.49328000    -1.07368900     0.00000000
O        1.17043500     1.16659000     0.00000000
H        2.97948800     0.25882900     0.00000000
H        0.49883300    -1.10719500     0.00000000
units angstrom
}
""", 0)

BENCH12_S22_HB_h2o_2 = input.process_input("""
molecule dimer {
0 1
O       -1.55100700    -0.11452000     0.00000000
H       -1.93425900     0.76250300     0.00000000
H       -0.59967700     0.04071200     0.00000000
O        1.35062500     0.11146900     0.00000000
H        1.68039800    -0.37374100    -0.75856100
H        1.68039800    -0.37374100     0.75856100
units angstrom
}
""", 0)

BENCH12_S22_HB_nh3_2 = input.process_input("""
molecule dimer {
0 1
N       -1.57871800    -0.04661100     0.00000000
H       -2.15862100     0.13639600    -0.80956500
H       -2.15862100     0.13639600     0.80956500
H       -0.84947100     0.65819300     0.00000000
N        1.57871800     0.04661100     0.00000000
H        2.15862100    -0.13639600    -0.80956500
H        0.84947100    -0.65819300     0.00000000
H        2.15862100    -0.13639600     0.80956500
units angstrom
}
""", 0)

BENCH12_S22_MX_bz_h2o = input.process_input("""
molecule dimer {
0 1
C        0.78061170    -0.60988750    -1.20754260
C        0.47840390     0.75104060    -1.20790400
C        0.32765920     1.43185730     0.00000000
C        0.47840390     0.75104060     1.20790400
C        0.78061170    -0.60988750     1.20754260
C        0.93215100    -1.28996140     0.00000000
H        0.89666880    -1.13760510    -2.14414820
H        0.35738950     1.27820910    -2.14405460
H        0.09185930     2.48714070     0.00000000
H        0.35738950     1.27820910     2.14405460
H        0.89666880    -1.13760510     2.14414820
H        1.16900640    -2.34516680     0.00000000
O       -2.78852700    -0.27448540     0.00000000
H       -2.62291140    -1.21908310     0.00000000
H       -1.90151030     0.09791100     0.00000000
units angstrom
}
""", 0)

BENCH12_S22_MX_bz_hcn = input.process_input("""
molecule dimer {
0 1
C       -0.70977410    -0.99042300     1.20770180
C       -1.40653400    -0.96535290     0.00000000
C       -0.70977410    -0.99042300    -1.20770180
C        0.68396510    -1.04051050    -1.20786520
C        1.38097790    -1.06555220     0.00000000
C        0.68396510    -1.04051050     1.20786520
H       -1.24994820    -0.96862800     2.14405070
H       -2.48691970    -0.92370600     0.00000000
H       -1.24994820    -0.96862800    -2.14405070
H        1.22428820    -1.05807530    -2.14425630
H        2.46158860    -1.10298180     0.00000000
H        1.22428820    -1.05807530     2.14425630
N       -0.00341180     3.53539260     0.00000000
C        0.07519630     2.37070400     0.00000000
H        0.14762950     1.30528470     0.00000000
units angstrom
}
""", 0)

BENCH12_S22_MX_bz_nh3 = input.process_input("""
molecule dimer {
0 1
C       -0.73928100     0.51587850    -1.20710790
C       -1.42614420     0.39654550     0.00000000
C       -0.73928100     0.51587850     1.20710790
C        0.63422690     0.75463980     1.20707350
C        1.32104340     0.87375660     0.00000000
C        0.63422690     0.75463980    -1.20707350
H       -1.27194950     0.42063160    -2.14328940
H       -2.49022050     0.20523810     0.00000000
H       -1.27194950     0.42063160     2.14328940
H        1.16680050     0.84748850     2.14369500
H        2.38635850     1.05963120     0.00000000
H        1.16680050     0.84748850    -2.14369500
N        0.18039300    -2.94912310     0.00000000
H        0.75954950    -3.14594770    -0.80607290
H        0.75954950    -3.14594770     0.80607290
H        0.04441670    -1.94493990     0.00000000
units angstrom
}
""", 0)

BENCH12_S22_MX_c2h2_c2h4 = input.process_input("""
molecule dimer {
0 1
C        0.00000000    -0.66757800    -2.12465900
C        0.00000000     0.66757800    -2.12465900
H        0.92362100    -1.23225300    -2.12618500
H       -0.92362100    -1.23225300    -2.12618500
H       -0.92362100     1.23225300    -2.12618500
H        0.92362100     1.23225300    -2.12618500
C        0.00000000     0.00000000     2.90050300
C        0.00000000     0.00000000     1.69324000
H        0.00000000     0.00000000     0.62735200
H        0.00000000     0.00000000     3.96392900
units angstrom
}
""", 0)

BENCH12_thio1 = input.process_input("""
molecule dimer {
0 1
S        0.00000000     0.00000000     0.00000000
C        0.00000000     0.00000000     1.72529600
H        0.93086973     0.00000000     2.26617608
C       -1.26880426    -0.00000000     2.22340435
H       -1.49163472     0.00003848     3.28002665
C       -2.25995860     0.00025938     1.20248869
H       -3.32270915     0.00043325     1.39397066
C       -1.72454170    -0.00024064    -0.05102823
H       -2.23765332    -0.00041661    -0.99748803
units angstrom
}
""", 0)

BENCH12_thio2 = input.process_input("""
molecule dimer {
0 1
S       -1.85234617    -1.20349129     0.00000000
C       -3.20012435    -0.12541479     0.00000000
H       -4.20210000    -0.51938625     0.00000000
C       -2.79844642     1.17619637     0.00000000
H       -3.48487528     2.00980042     0.00000000
C       -1.38764366     1.31782604     0.00000000
H       -0.88527164     2.27414792     0.00000000
C       -0.71318040     0.12170249     0.00000000
C        0.71318040    -0.12170249     0.00000000
C        1.38764366    -1.31782604     0.00000000
H        0.88527164    -2.27414792     0.00000000
C        2.79844642    -1.17619637     0.00000000
H        3.48487528    -2.00980042     0.00000000
C        3.20012435     0.12541479     0.00000000
H        4.20210000     0.51938625     0.00000000
S        1.85234617     1.20349129     0.00000000
units angstrom
}
""", 0)

# <<< Geometry Specification Strings >>>
rxnpattern = re.compile(r'^(.+)-(.+)-(.+)$')
GEOS = {}
for rxn in HRXN:
    for rgt in ACTV['%s-%s' % (dbse, rxn)]:

        molname = rxnpattern.match(rgt)
        GEOS['%s' % (rgt)] = eval('%s_%s' % (dbse, molname.group(2)))
