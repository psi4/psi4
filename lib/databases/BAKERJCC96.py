import re
import input

# <<< BAKERJCC96 Database Module >>>
# Geometries from Baker and Chan J. Comput. Chem. 17 888 (1996), as reported
# in Bakken and Helgaker, J. Chem. Phys. 117, 9160 (2002).
dbse = 'BAKERJCC96'
isOS = 'true'

# <<< Database Members >>>
HRXN = ['HCN_to_HNC', 'HCCH_to_CCH2', 'H2CO_to_H2_CO', 'parent_diels_alder', 's_tetrazine_to_2HCN_N2', 'CH3CH3_to_CH2CH2_H2', 'CH3CH2F_to_CH2CH2_HF', 'CH2CHOH_to_CH3CHO', 'silylene_insertion', 'HNCCS_to_HCN_CS', 'acrolein_rotation', 'HCONHOH_to_HCOHNHO', 'HNC_H2_to_H2CNH', 'H2CNH_to_HCNH2', 'HCNH2_to_HCN_H2']
HRXN_TEMP = ['parent_diels_alder']

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV['%s-%s'            % (dbse, 'HCN_to_HNC' )] = ['%s-%s-reagent'      % (dbse, 'HCN_to_HNC')]
RXNM['%s-%s'            % (dbse, 'HCN_to_HNC' )] = dict(zip(ACTV['%s-%s' % (dbse, 'HCN_to_HNC')], [+1]))

ACTV['%s-%s'            % (dbse, 'HCCH_to_CCH2' )] = ['%s-%s-reagent'      % (dbse, 'HCCH_to_CCH2')]
RXNM['%s-%s'            % (dbse, 'HCCH_to_CCH2' )] = dict(zip(ACTV['%s-%s' % (dbse, 'HCCH_to_CCH2')], [+1]))

ACTV['%s-%s'            % (dbse, 'H2CO_to_H2_CO' )] = ['%s-%s-reagent'      % (dbse, 'H2CO_to_H2_CO')]
RXNM['%s-%s'            % (dbse, 'H2CO_to_H2_CO' )] = dict(zip(ACTV['%s-%s' % (dbse, 'H2CO_to_H2_CO')], [+1]))

ACTV['%s-%s'            % (dbse, 'parent_diels_alder' )] = ['%s-%s-reagent'      % (dbse, 'parent_diels_alder')]
RXNM['%s-%s'            % (dbse, 'parent_diels_alder' )] = dict(zip(ACTV['%s-%s' % (dbse, 'parent_diels_alder')], [+1]))

ACTV['%s-%s'            % (dbse, 's_tetrazine_to_2HCN_N2' )] = ['%s-%s-reagent'      % (dbse, 's_tetrazine_to_2HCN_N2')]
RXNM['%s-%s'            % (dbse, 's_tetrazine_to_2HCN_N2' )] = dict(zip(ACTV['%s-%s' % (dbse, 's_tetrazine_to_2HCN_N2')], [+1]))

ACTV['%s-%s'            % (dbse, 'CH3CH3_to_CH2CH2_H2' )] = ['%s-%s-reagent'      % (dbse, 'CH3CH3_to_CH2CH2_H2')]
RXNM['%s-%s'            % (dbse, 'CH3CH3_to_CH2CH2_H2' )] = dict(zip(ACTV['%s-%s' % (dbse, 'CH3CH3_to_CH2CH2_H2')], [+1]))

ACTV['%s-%s'            % (dbse, 'CH3CH2F_to_CH2CH2_HF' )] = ['%s-%s-reagent'      % (dbse, 'CH3CH2F_to_CH2CH2_HF')]
RXNM['%s-%s'            % (dbse, 'CH3CH2F_to_CH2CH2_HF' )] = dict(zip(ACTV['%s-%s' % (dbse, 'CH3CH2F_to_CH2CH2_HF')], [+1]))

ACTV['%s-%s'            % (dbse, 'CH2CHOH_to_CH3CHO' )] = ['%s-%s-reagent'      % (dbse, 'CH2CHOH_to_CH3CHO')]
RXNM['%s-%s'            % (dbse, 'CH2CHOH_to_CH3CHO' )] = dict(zip(ACTV['%s-%s' % (dbse, 'CH2CHOH_to_CH3CHO')], [+1]))

ACTV['%s-%s'            % (dbse, 'silylene_insertion' )] = ['%s-%s-reagent'      % (dbse, 'silylene_insertion')]
RXNM['%s-%s'            % (dbse, 'silylene_insertion' )] = dict(zip(ACTV['%s-%s' % (dbse, 'silylene_insertion')], [+1]))

ACTV['%s-%s'            % (dbse, 'HNCCS_to_HCN_CS' )] = ['%s-%s-reagent'      % (dbse, 'HNCCS_to_HCN_CS')]
RXNM['%s-%s'            % (dbse, 'HNCCS_to_HCN_CS' )] = dict(zip(ACTV['%s-%s' % (dbse, 'HNCCS_to_HCN_CS')], [+1]))

ACTV['%s-%s'            % (dbse, 'acrolein_rotation' )] = ['%s-%s-reagent'      % (dbse, 'acrolein_rotation')]
RXNM['%s-%s'            % (dbse, 'acrolein_rotation' )] = dict(zip(ACTV['%s-%s' % (dbse, 'acrolein_rotation')], [+1]))

ACTV['%s-%s'            % (dbse, 'HCONHOH_to_HCOHNHO' )] = ['%s-%s-reagent'      % (dbse, 'HCONHOH_to_HCOHNHO')]
RXNM['%s-%s'            % (dbse, 'HCONHOH_to_HCOHNHO' )] = dict(zip(ACTV['%s-%s' % (dbse, 'HCONHOH_to_HCOHNHO')], [+1]))

ACTV['%s-%s'            % (dbse, 'HNC_H2_to_H2CNH' )] = ['%s-%s-reagent'      % (dbse, 'HNC_H2_to_H2CNH')]
RXNM['%s-%s'            % (dbse, 'HNC_H2_to_H2CNH' )] = dict(zip(ACTV['%s-%s' % (dbse, 'HNC_H2_to_H2CNH')], [+1]))

ACTV['%s-%s'            % (dbse, 'H2CNH_to_HCNH2' )] = ['%s-%s-reagent'      % (dbse, 'H2CNH_to_HCNH2')]
RXNM['%s-%s'            % (dbse, 'H2CNH_to_HCNH2' )] = dict(zip(ACTV['%s-%s' % (dbse, 'H2CNH_to_HCNH2')], [+1]))

ACTV['%s-%s'            % (dbse, 'HCNH2_to_HCN_H2' )] = ['%s-%s-reagent'      % (dbse, 'HCNH2_to_HCN_H2')]
RXNM['%s-%s'            % (dbse, 'HCNH2_to_HCN_H2' )] = dict(zip(ACTV['%s-%s' % (dbse, 'HCNH2_to_HCN_H2')], [+1]))

# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'            % (dbse, 'HCN_to_HNC'             )] =    0.000
BIND['%s-%s'            % (dbse, 'HCCH_to_CCH2'           )] =    0.000
BIND['%s-%s'            % (dbse, 'H2CO_to_H2_CO'          )] =    0.000
BIND['%s-%s'            % (dbse, 'parent_diels_alder'     )] =    0.000
BIND['%s-%s'            % (dbse, 's_tetrazine_to_2HCN_N2' )] =    0.000
BIND['%s-%s'            % (dbse, 'CH3CH3_to_CH2CH2_H2'    )] =    0.000
BIND['%s-%s'            % (dbse, 'CH3CH2F_to_CH2CH2_HF'   )] =    0.000
BIND['%s-%s'            % (dbse, 'CH2CHOH_to_CH3CHO'      )] =    0.000
BIND['%s-%s'            % (dbse, 'silylene_insertion'     )] =    0.000
BIND['%s-%s'            % (dbse, 'HNCCS_to_HCN_CS'        )] =    0.000
BIND['%s-%s'            % (dbse, 'acrolein_rotation'      )] =    0.000
BIND['%s-%s'            % (dbse, 'HCONHOH_to_HCOHNHO'     )] =    0.000
BIND['%s-%s'            % (dbse, 'HNC_H2_to_H2CNH'        )] =    0.000
BIND['%s-%s'            % (dbse, 'H2CNH_to_HCNH2'         )] =    0.000
BIND['%s-%s'            % (dbse, 'HCNH2_to_HCN_H2'        )] =    0.000

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, 'HCN_to_HNC' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'HCN_to_HNC' )] = ''
TAGL['%s-%s'            % (dbse, 'HCCH_to_CCH2' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'HCCH_to_CCH2' )] = ''
TAGL['%s-%s'            % (dbse, 'H2CO_to_H2_CO' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'H2CO_to_H2_CO' )] = ''
TAGL['%s-%s'            % (dbse, 'parent_diels_alder' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'parent_diels_alder' )] = ''
TAGL['%s-%s'            % (dbse, 's_tetrazine_to_2HCN_N2' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 's_tetrazine_to_2HCN_N2' )] = ''
TAGL['%s-%s'            % (dbse, 'CH3CH3_to_CH2CH2_H2' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'CH3CH3_to_CH2CH2_H2' )] = ''
TAGL['%s-%s'            % (dbse, 'CH3CH2F_to_CH2CH2_HF' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'CH3CH2F_to_CH2CH2_HF' )] = ''
TAGL['%s-%s'            % (dbse, 'CH2CHOH_to_CH3CHO' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'CH2CHOH_to_CH3CHO' )] = ''
TAGL['%s-%s'            % (dbse, 'silylene_insertion' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'silylene_insertion' )] = ''
TAGL['%s-%s'            % (dbse, 'HNCCS_to_HCN_CS' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'HNCCS_to_HCN_CS' )] = ''
TAGL['%s-%s'            % (dbse, 'acrolein_rotation' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'acrolein_rotation' )] = ''
TAGL['%s-%s'            % (dbse, 'HCONHOH_to_HCOHNHO' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'HCONHOH_to_HCOHNHO' )] = ''
TAGL['%s-%s'            % (dbse, 'HNC_H2_to_H2CNH' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'HNC_H2_to_H2CNH' )] = ''
TAGL['%s-%s'            % (dbse, 'H2CNH_to_HCNH2' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'H2CNH_to_HCNH2' )] = ''
TAGL['%s-%s'            % (dbse, 'HCNH2_to_HCN_H2' )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'HCNH2_to_HCN_H2' )] = ''

# <<< Molecule Specifications >>>
BAKERJCC96_HCN_to_HNC = input.process_input("""
molecule {
0 1
C  -0.0399606537   1.7574844925   0.0000000000
N   1.2331003808   0.0000000000   0.0000000000
H  -1.1931397271  -1.7574844925   0.0000000000
units bohr
}
set { guess gwh }
""")

BAKERJCC96_HCCH_to_CCH2 = input.process_input("""
molecule {
0 1
C  -0.4287449922  -0.0396754553   0.0000000000
C   1.9151999208   0.0000000000   0.0000000000
H  -2.1914826494   0.9243341981   0.0000000000
H   1.1657551200  -2.7344386158   0.0000000000
units bohr
}
""")

BAKERJCC96_H2CO_to_H2_CO = input.process_input("""
molecule {
0 1
C   0.4656871259   0.8485069310   0.0000000000
O   2.6701879709   0.0000000000   0.0000000000
H  -0.8014735829  -1.2561076240   0.0000000000
H  -2.3344015139   0.4076006930   0.0000000000
units bohr
}
""")

BAKERJCC96_parent_diels_alder = input.process_input("""
molecule {
0 1
C   2.8332856188   0.0000000000   1.3228081920
C   2.8332856188   0.0000000000  -1.3228081920
C  -2.3253041766   1.1773255830   1.4550890112
C  -2.3253041766   1.1773255830  -1.4550890112
C  -0.9003485915  -0.5867135859   2.6456163840
C  -0.9003485915  -0.5867135859  -2.6456163840
H  -3.5271811018   2.4200543219   2.5397117805
H  -3.5271811018   2.4200543219  -2.5397117805
H   3.2953993811  -1.7330277495   2.2966434467
H   3.2953993811  -1.7330277495  -2.2966434467
H   2.7098759932   1.7961080923   2.2800674415
H   2.7098759932   1.7961080923  -2.2800674415
H   0.3094408240  -1.7857534614   1.5179110011
H   0.3094408240  -1.7857534614  -1.5179110011
H  -1.0137937565  -0.6962109120   4.6804210291
H  -1.0137937565  -0.6962109120  -4.6804210291
units bohr
}
""")

BAKERJCC96_s_tetrazine_to_2HCN_N2 = input.process_input("""
molecule {
0 1
N   1.4172944914   2.0866021582   0.0000000000
N  -1.4172944914   2.0866021582   0.0000000000
N   1.1338355931  -2.3320847650   0.0000000000
N  -1.1338355931  -2.3320847650   0.0000000000
C   2.5511300846   0.1227413034   0.0000000000
C  -2.5511300846   0.1227413034   0.0000000000
H   4.5920341522   0.1227413034   0.0000000000
H  -4.5920341522   0.1227413034   0.0000000000
units bohr
}
""")

BAKERJCC96_CH3CH3_to_CH2CH2_H2 = input.process_input("""
molecule {
0 1
C  -0.8859807599  -0.5993092127   0.0000000000
C   1.7490330691   0.0000000000   0.0000000000
H  -2.8124044495   1.4800624400   0.0000000000
H   2.6980243606  -0.2473177820   1.8113616294
H   2.6980243606  -0.2473177820  -1.8113616294
H  -1.6346529801  -1.2327446580   1.8113616294
H  -1.6346529801  -1.2327446580  -1.8113616294
H  -0.1773906205   2.0793716527   0.0000000000
units bohr
}
""")

BAKERJCC96_CH3CH2F_to_CH2CH2_HF = input.process_input("""
molecule {
0 1
C  -1.2114041996   0.8408020855   0.0000000000
C   0.2710949367  -1.4185486447   0.0000000000
F   3.5694675519   0.0000000000   0.0000000000
H   0.8525043764   2.7837833236   0.0000000000
H   0.3440183107  -2.4023246321   1.7866606845
H   0.3440183107  -2.4023246321  -1.7866606845
H  -2.0848496434   1.2993062499   1.7866606845
H  -2.0848496434   1.2993062499  -1.7866606845
units bohr
}
""")

BAKERJCC96_CH2CHOH_to_CH3CHO = input.process_input("""
molecule {
0 1
C  -0.8638822546   0.0813052799  -1.2005986026
C   0.8568384546   0.0000000000   0.8814626606
O   0.0000000000   0.0000000000   3.1838367722
H  -0.6251819034   1.8414643222  -2.2435002266
H  -1.6228517906  -0.1159760771   1.2188517026
H   2.8765713088  -0.2172341086   0.5404787182
H  -0.6214938148  -1.5895594164  -2.3805310244
units bohr
}
""")

BAKERJCC96_silylene_insertion = input.process_input("""
molecule {
0 1
C  -1.0718066965  -0.1928298581  -2.5299430030
C   1.0250683554   0.0000000000  -0.5763517845
Si  0.0000000000   0.0000000000   3.6474101733
H  -2.1062850458  -1.9954618846  -2.4930629368
H   2.3105025336   0.4726358860   5.1508238324
H  -0.1255963623  -1.7429478038   1.4155380424
H  -1.8281836996   2.0590612062   4.1373801504
H   2.1120964044   1.7375623482  -0.9231008641
H   2.3199575761  -1.5749557960  -0.8688541874
H  -0.2137259723  -0.1257146643  -4.4428207722
H  -2.4220270929   1.3626505664  -2.5170186506
units bohr
}
""")

BAKERJCC96_HNCCS_to_HCN_CS = input.process_input("""
molecule {
0 1
H  -3.8046808257   0.6512707943   0.0000000000
N  -2.5162100514  -0.7568079564   0.0000000000
C  -0.3115625330  -1.2877146607   0.0000000000
C   2.4597092296   0.8743461413   0.0000000000
S   5.2948168550   0.0000000000   0.0000000000
units bohr
}
""")

BAKERJCC96_acrolein_rotation = input.process_input("""
molecule {
0 1
C  -1.2039490098   0.6299819763  -1.8443710751
C   1.2056120532   0.0000000000   1.9839551839
C   0.0666798969  -0.9449729645  -0.3221285864
O   0.0000000000   0.0000000000   3.9490688445
H  -2.0522570986  -0.0738599558  -3.5620058833
H   3.1211832712   0.7038419322   1.9619811778
H   0.2424607776  -2.9181815653  -0.8128813887
H  -1.3797298905   2.6031905771  -1.3536182728
units bohr
}
""")

BAKERJCC96_HCONHOH_to_HCOHNHO = input.process_input("""
molecule {
0 1
O  -3.5099961475   0.0808581274   0.0000000000
O   3.0085718011   0.0000000000   0.0000000000
C  -1.2730134573   1.0962433404   0.0000000000
N   0.4462296750  -0.6585505903   0.0000000000
H  -0.7144888327   3.0985015808   0.0000000000
H  -2.0229706827  -2.0506512675   0.0000000000
H   4.0656676440  -1.5664011909   0.0000000000
units bohr
}
""")

BAKERJCC96_HNC_H2_to_H2CNH = input.process_input("""
molecule {
0 1
H   1.3134432637  -0.2040854069   3.0062409561
H  -2.0477645147   0.2508396893  -1.7949671795
H   0.1303083495  -0.0467542824  -2.3515361651
N   0.0000000000   0.0000000000   1.6630059611
C   0.6040129016   0.0000000000  -0.5227435726
units bohr
}
""")

BAKERJCC96_H2CNH_to_HCNH2 = input.process_input("""
molecule {
0 1
H  -1.5502016034   1.0481518658   1.5179225979
H   0.3563149954   1.0481518658  -2.1765524295
H   0.1539684632  -2.0963037316   0.1645151437
N   0.0000000000   0.0000000000   1.2546414423
C   1.0399181448   0.0000000000  -0.7605267545
units bohr
}
""")

BAKERJCC96_HCNH2_to_HCN_H2 = input.process_input("""
molecule {
0 1
C   1.4474353774   0.0000000000  -1.1288243208
N   0.0000000000   0.0000000000   0.9719363835
H   3.2098577258   0.2829784675  -0.5084575223
H  -2.3053050298  -1.2021721464   0.4824091612
H  -2.3519880735   0.9191936788   0.1829362984
units bohr
}
""")


# <<< Geometry Specification Strings >>>
rxnpattern = re.compile(r'^(.+)-(.+)-(.+)$')
GEOS = {}
for rxn in HRXN:
   for rgt in ACTV['%s-%s' % (dbse, rxn)]:

            molname = rxnpattern.match(rgt)
            GEOS['%s' % (rgt)] = eval('%s_%s' % (dbse, molname.group(2)))

