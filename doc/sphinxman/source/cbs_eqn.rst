.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2017 The Psi4 Developers.
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

.. math:: E_{\text{total}}^{\text{CBS}} = \mathcal{F}_{\textbf{scf_scheme}} \left(E_{\text{total},\; \text{SCF}}^{\textbf{scf_basis}}\right) \; + \mathcal{F}_{\textbf{corl_scheme}} \left(E_{\text{corl},\; \textbf{corl_wfn}}^{\textbf{corl_basis}}\right) \; + \delta_{\textbf{delta_wfn_lesser}}^{\textbf{delta_wfn}} \; + \delta_{\textbf{delta2_wfn_lesser}}^{\textbf{delta2_wfn}} \; + \delta_{\textbf{delta3_wfn_lesser}}^{\textbf{delta3_wfn}} \; + \delta_{\textbf{delta4_wfn_lesser}}^{\textbf{delta4_wfn}} \; + \delta_{\textbf{delta5_wfn_lesser}}^{\textbf{delta5_wfn}}

Here, :math:`\mathcal{F}` is an energy or energy extrapolation scheme, and the following also hold.

.. math:: \delta_{\textbf{delta_wfn_lesser}}^{\textbf{delta_wfn}} \; = \mathcal{F}_{\textbf{delta_scheme}} \left(E_{\text{corl},\; \textbf{delta_wfn}}^{\textbf{delta_basis}}\right) - \mathcal{F}_{\textbf{delta_scheme}} \left(E_{\text{corl},\; \textbf{delta_wfn_lesser}}^{\textbf{delta_basis}}\right)

.. math:: \delta_{\textbf{delta2_wfn_lesser}}^{\textbf{delta2_wfn}} \; = \mathcal{F}_{\textbf{delta2_scheme}} \left(E_{\text{corl},\; \textbf{delta2_wfn}}^{\textbf{delta2_basis}}\right) - \mathcal{F}_{\textbf{delta2_scheme}} \left(E_{\text{corl},\; \textbf{delta2_wfn_lesser}}^{\textbf{delta2_basis}}\right)

.. math:: \delta_{\textbf{delta3_wfn_lesser}}^{\textbf{delta3_wfn}} \; = \mathcal{F}_{\textbf{delta3_scheme}} \left(E_{\text{corl},\; \textbf{delta3_wfn}}^{\textbf{delta3_basis}}\right) - \mathcal{F}_{\textbf{delta3_scheme}} \left(E_{\text{corl},\; \textbf{delta3_wfn_lesser}}^{\textbf{delta3_basis}}\right)

.. math:: \delta_{\textbf{delta4_wfn_lesser}}^{\textbf{delta4_wfn}} \; = \mathcal{F}_{\textbf{delta4_scheme}} \left(E_{\text{corl},\; \textbf{delta4_wfn}}^{\textbf{delta4_basis}}\right) - \mathcal{F}_{\textbf{delta4_scheme}} \left(E_{\text{corl},\; \textbf{delta4_wfn_lesser}}^{\textbf{delta4_basis}}\right)

.. math:: \delta_{\textbf{delta5_wfn_lesser}}^{\textbf{delta5_wfn}} \; = \mathcal{F}_{\textbf{delta5_scheme}} \left(E_{\text{corl},\; \textbf{delta5_wfn}}^{\textbf{delta5_basis}}\right) - \mathcal{F}_{\textbf{delta5_scheme}} \left(E_{\text{corl},\; \textbf{delta5_wfn_lesser}}^{\textbf{delta5_basis}}\right)

A translation of this ungainly equation to example [5] below is as
follows. In words, this is a double- and triple-zeta 2-point
Helgaker-extrapolated CCSD(T) coupled-cluster correlation correction
appended to a triple- and quadruple-zeta 2-point
Helgaker-extrapolated MP2 correlation energy appended to a SCF/aug-cc-pVQZ
reference energy.

.. math:: E_{\text{total}}^{\text{CBS}} = \mathcal{F}_{\text{highest_1}} \left(E_{\text{total},\; \text{SCF}}^{\text{aug-cc-pVQZ}}\right) \; + \mathcal{F}_{\text{corl_xtpl_helgaker_2}} \left(E_{\text{corl},\; \text{MP2}}^{\text{aug-cc-pV[TQ]Z}}\right) \; + \delta_{\text{MP2}}^{\text{CCSD(T)}}

.. math:: \delta_{\text{MP2}}^{\text{CCSD(T)}} \; = \mathcal{F}_{\text{corl_xtpl_helgaker_2}} \left(E_{\text{corl},\; \text{CCSD(T)}}^{\text{aug-cc-pV[DT]Z}}\right) - \mathcal{F}_{\text{corl_xtpl_helgaker_2}} \left(E_{\text{corl},\; \text{MP2}}^{\text{aug-cc-pV[DT]Z}}\right)


