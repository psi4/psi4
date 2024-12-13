# Objectives
- [ ] Add external potential to SAPT(DFT):
- [ ] Derive FSAPT for SAPT(DFT) equations. Maybe there are algorithm choices
        for FSAPT implementation that are actually faster?

# SAPT(DFT) - External Potential
1. use scf_helper() for handling external potential
    - shouldn't need to really do anything here now that initial_cartesians_
    are set?
2. Need to extract out external potential from FISAPT::nuclear():
    1. `/home/amwalla3/gits/psi4/psi4/driver/procrouting/sapt/sapt_jk_terms.py:75`  
       add result to V_A and V_B and should propagate throughout SAPT(DFT)
```python
    # psi4/driver/procrouting/sapt/sapt_jk_terms.py
    # Potential ints
    mints = core.MintsHelper(wfn_A.basisset())
    cache["V_A"] = mints.ao_potential()
    # TODO: add external potential
    # cache["V_A"].axpy(1.0, wfn_A.Va())

    mints = core.MintsHelper(wfn_B.basisset())
    cache["V_B"] = mints.ao_potential()
```
3. Electrostatics energy:
    1. Add something similar to:
```cpp
// psi4/src/psi4/fisapt/fisapt.cc
// External potential interactions
if (reference_->has_potential_variable("A") && reference_->has_potential_variable("B")) {
// Add the interaction between the external potentials in A and B
// Multiply by 2 to get the full A-B + B-A interaction energy
scalars_["Extern-Extern"] = matrices_["extern_extern_IE"]->get(0, 1)*2.0; 
outfile->Printf("    Extern-Extern       = %18.12lf [Eh]\n", scalars_["Extern-Extern"]);
}
```
4. felst() can wait for now in SAPT(DFT) until FSAPT(DFT) formulation figured
   out...
