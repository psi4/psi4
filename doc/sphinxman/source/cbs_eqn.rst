
.. index:: 
   pair: cbs(); theory

.. math:: E_{total}^{\text{CBS}} = \mathcal{F}_{\textbf{scf\_scheme}} \left(E_{total,\; \text{SCF}}^{\textbf{scf\_basis}}\right) \; + \mathcal{F}_{\textbf{corl\_scheme}} \left(E_{corl,\; \textbf{corl\_wfn}}^{\textbf{corl\_basis}}\right) \; + \delta_{\textbf{delta\_wfn\_lesser}}^{\textbf{delta\_wfn}} \; + \delta_{\textbf{delta2\_wfn\_lesser}}^{\textbf{delta2\_wfn}}

Here, :math:`\mathcal{F}` is an energy or energy extrapolation scheme, and the following also hold.

.. math:: \delta_{\textbf{delta\_wfn\_lesser}}^{\textbf{delta\_wfn}} \; = \mathcal{F}_{\textbf{delta\_scheme}} \left(E_{corl,\; \textbf{delta\_wfn}}^{\textbf{delta\_basis}}\right) - \mathcal{F}_{\textbf{delta\_scheme}} \left(E_{corl,\; \textbf{delta\_wfn\_lesser}}^{\textbf{delta\_basis}}\right)

.. math:: \delta_{\textbf{delta2\_wfn\_lesser}}^{\textbf{delta2\_wfn}} \; = \mathcal{F}_{\textbf{delta2\_scheme}} \left(E_{corl,\; \textbf{delta2\_wfn}}^{\textbf{delta2\_basis}}\right) - \mathcal{F}_{\textbf{delta2\_scheme}} \left(E_{corl,\; \textbf{delta2\_wfn\_lesser}}^{\textbf{delta2\_basis}}\right)

A translation of this ungainly equation to example [5] below is as
follows. In words, this is a double- and triple-zeta 2-point
Helgaker-extrapolated CCSD(T) coupled-cluster correlation correction
appended to a triple- and quadruple-zeta 2-point
Helgaker-extrapolated MP2 correlation energy appended to a SCF/aug-cc-pVQZ
reference energy.

.. math:: E_{total}^{\text{CBS}} = \mathcal{F}_{\text{highest\_1}} \left(E_{total,\; \text{SCF}}^{\text{aug-cc-pVQZ}}\right) \; + \mathcal{F}_{\text{corl\_xtpl\_helgaker\_2}} \left(E_{corl,\; \text{MP2}}^{\text{aug-cc-pV[TQ]Z}}\right) \; + \delta_{\text{MP2}}^{\text{CCSD(T)}}

.. math:: \delta_{\text{MP2}}^{\text{CCSD(T)}} \; = \mathcal{F}_{\text{corl\_xtpl\_helgaker\_2}} \left(E_{corl,\; \text{CCSD(T)}}^{\text{aug-cc-pV[DT]Z}}\right) - \mathcal{F}_{\text{corl\_xtpl\_helgaker\_2}} \left(E_{corl,\; \text{MP2}}^{\text{aug-cc-pV[DT]Z}}\right)


