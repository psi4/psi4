.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2016 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This program is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU General Public License as published by
.. # the Free Software Foundation; either version 2 of the License, or
.. # (at your option) any later version.
.. #
.. # This program is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU General Public License for more details.
.. #
.. # You should have received a copy of the GNU General Public License along
.. # with this program; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #

.. index:: 
   pair: cbs(); theory

.. _`eq:cbs`:

.. math:: E_{total}^{\text{CBS}} = \mathcal{F}_{\textbf{scf\_scheme}} \left(E_{total,\; \text{SCF}}^{\textbf{scf\_basis}}\right) \; + \mathcal{F}_{\textbf{corl\_scheme}} \left(E_{corl,\; \textbf{corl\_wfn}}^{\textbf{corl\_basis}}\right) \; + \delta_{\textbf{delta\_wfn\_lesser}}^{\textbf{delta\_wfn}} \; + \delta_{\textbf{delta2\_wfn\_lesser}}^{\textbf{delta2\_wfn}} \; + \delta_{\textbf{delta3\_wfn\_lesser}}^{\textbf{delta3\_wfn}} \; + \delta_{\textbf{delta4\_wfn\_lesser}}^{\textbf{delta4\_wfn}} \; + \delta_{\textbf{delta5\_wfn\_lesser}}^{\textbf{delta5\_wfn}}

Here, :math:`\mathcal{F}` is an energy or energy extrapolation scheme, and the following also hold.

.. math:: \delta_{\textbf{delta\_wfn\_lesser}}^{\textbf{delta\_wfn}} \; = \mathcal{F}_{\textbf{delta\_scheme}} \left(E_{corl,\; \textbf{delta\_wfn}}^{\textbf{delta\_basis}}\right) - \mathcal{F}_{\textbf{delta\_scheme}} \left(E_{corl,\; \textbf{delta\_wfn\_lesser}}^{\textbf{delta\_basis}}\right)

.. math:: \delta_{\textbf{delta2\_wfn\_lesser}}^{\textbf{delta2\_wfn}} \; = \mathcal{F}_{\textbf{delta2\_scheme}} \left(E_{corl,\; \textbf{delta2\_wfn}}^{\textbf{delta2\_basis}}\right) - \mathcal{F}_{\textbf{delta2\_scheme}} \left(E_{corl,\; \textbf{delta2\_wfn\_lesser}}^{\textbf{delta2\_basis}}\right)

.. math:: \delta_{\textbf{delta3\_wfn\_lesser}}^{\textbf{delta3\_wfn}} \; = \mathcal{F}_{\textbf{delta3\_scheme}} \left(E_{corl,\; \textbf{delta3\_wfn}}^{\textbf{delta3\_basis}}\right) - \mathcal{F}_{\textbf{delta3\_scheme}} \left(E_{corl,\; \textbf{delta3\_wfn\_lesser}}^{\textbf{delta3\_basis}}\right)

.. math:: \delta_{\textbf{delta4\_wfn\_lesser}}^{\textbf{delta4\_wfn}} \; = \mathcal{F}_{\textbf{delta4\_scheme}} \left(E_{corl,\; \textbf{delta4\_wfn}}^{\textbf{delta4\_basis}}\right) - \mathcal{F}_{\textbf{delta4\_scheme}} \left(E_{corl,\; \textbf{delta4\_wfn\_lesser}}^{\textbf{delta4\_basis}}\right)

.. math:: \delta_{\textbf{delta5\_wfn\_lesser}}^{\textbf{delta5\_wfn}} \; = \mathcal{F}_{\textbf{delta5\_scheme}} \left(E_{corl,\; \textbf{delta5\_wfn}}^{\textbf{delta5\_basis}}\right) - \mathcal{F}_{\textbf{delta5\_scheme}} \left(E_{corl,\; \textbf{delta5\_wfn\_lesser}}^{\textbf{delta5\_basis}}\right)

A translation of this ungainly equation to example [5] below is as
follows. In words, this is a double- and triple-zeta 2-point
Helgaker-extrapolated CCSD(T) coupled-cluster correlation correction
appended to a triple- and quadruple-zeta 2-point
Helgaker-extrapolated MP2 correlation energy appended to a SCF/aug-cc-pVQZ
reference energy.

.. math:: E_{total}^{\text{CBS}} = \mathcal{F}_{\text{highest\_1}} \left(E_{total,\; \text{SCF}}^{\text{aug-cc-pVQZ}}\right) \; + \mathcal{F}_{\text{corl\_xtpl\_helgaker\_2}} \left(E_{corl,\; \text{MP2}}^{\text{aug-cc-pV[TQ]Z}}\right) \; + \delta_{\text{MP2}}^{\text{CCSD(T)}}

.. math:: \delta_{\text{MP2}}^{\text{CCSD(T)}} \; = \mathcal{F}_{\text{corl\_xtpl\_helgaker\_2}} \left(E_{corl,\; \text{CCSD(T)}}^{\text{aug-cc-pV[DT]Z}}\right) - \mathcal{F}_{\text{corl\_xtpl\_helgaker\_2}} \left(E_{corl,\; \text{MP2}}^{\text{aug-cc-pV[DT]Z}}\right)


