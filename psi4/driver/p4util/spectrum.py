#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

__all__ = [
    "Gaussian",
    "Lineshape",
    "Lorentzian",
    "prefactor_ecd",
    "prefactor_opa",
    "spectrum",
]

from abc import abstractmethod
from dataclasses import dataclass
from typing import Callable, Dict, List, Union

import numpy as np

from ..constants import constants


@dataclass
class Lineshape:
    """Lineshape ABC

    Notes
    -----
    domain -- Domain of the spectral band.
    gamma -- A function returning the broadening factor.
    Why do we use a callable broadening factor?
    For plots in the *wavelength domain*, the broadening factor depends on the location of the band's maximum.
    """
    domain: Union[np.ndarray, List[float]]
    gamma: Callable[[float], float]

    @abstractmethod
    def lineshape(self, x_0: float) -> np.ndarray:
        pass

    @abstractmethod
    def maximum(self, x_0: float) -> float:
        pass


class Gaussian(Lineshape):
    r"""Gaussian function on `domain`, centered at `x_0` with broadening `gamma`.

    Parameters
    ----------
    domain
        The domain of the Gaussian profile.
    gamma
        Broadening parameter.
        This is related to the full width at half maximum as :math:`\mathrm{FWHM} = \gamma \sqrt{2\ln 2}`

    Notes
    -----
    Use this profile to model inhomegenous broadening.
    """

    def lineshape(self, x_0: float) -> np.ndarray:
        """Gaussian function on :py:attr:`Lineshape.domain`, centered at `x_0` with broadening :py:attr:`Lineshape.gamma`.

        Parameters
        ----------
        x_0
            Center of the Gaussian, i.e. its maximum.

        Returns
        -------
        numpy.ndarray
            The Gaussian profile.

        """
        prefactor = 2.0 / (self.gamma(x_0) * np.sqrt(2.0 * np.pi))
        exponent = -2.0 * ((self.domain - x_0) / self.gamma(x_0))**2

        return prefactor * np.exp(exponent)

    def maximum(self, x_0: float) -> float:
        """Maximum value of Gaussian profile centered at `x_0`.

        Parameters
        ----------
        x_0
            Center of the Lorentzian, i.e. its maximum.

        """
        return 2.0 / (self.gamma(x_0) * np.sqrt(2.0 * np.pi))


class Lorentzian(Lineshape):
    """Lorentzian function on `domain`, centered at `x_0` with broadening `gamma`.

    Parameters
    ----------
    domain
        The domain of the Lorentzian profile.
    gamma
        Broadening parameter.
        This is the full width at half maximum (FWHM).

    Notes
    -----
    Use this profile to model homogeneous broadening.
    """

    def lineshape(self, x_0: float) -> np.ndarray:
        """Lorentzian function on :py:attr:`Lineshape.domain`, centered at `x_0` with broadening :py:attr:`Lineshape.gamma`.

        Parameters
        ----------
        x_0
            Center of the Lorentzian, i.e. its maximum.

        Returns
        -------
        numpy.ndarray
            The Lorentzian profile.
        """
        prefactor = 1.0 / np.pi
        numerator = self.gamma(x_0) / 2.0
        denominator = (self.domain - x_0)**2 + numerator**2

        return prefactor * (numerator / denominator)

    def maximum(self, x_0: float) -> float:
        """Maximum value of Lorentzian profile centered at `x_0`.

        Parameters
        ----------
        x_0
            Center of the Lorentzian, i.e. its maximum.

        """
        return 2.0 / (np.pi * self.gamma(x_0))


def prefactor_opa() -> float:
    r"""Prefactor for converting microscopic observable to decadic molar
    extinction coefficient in one-photon absorption.

    Notes
    -----
    This function implements the calculation of the following prefactor:

    .. math::

        k = \frac{4\pi^{2}N_{\mathrm{A}}}{3\times 1000\times \ln(10) (4 \pi \epsilon_{0}) n \hbar c}

    The prefactor is computed in SI units and then adjusted for the fact that
    we use atomic units to express microscopic observables: excitation energies
    and transition dipole moments.
    The refractive index :math:`n` is, in general, frequency-dependent. We
    assume it to be constant and equal to 1.
    """

    N_A = constants.get("Avogadro constant")
    c = constants.get("speed of light in vacuum")
    hbar = constants.get("Planck constant over 2 pi")
    e_0 = constants.get("electric constant")
    au_to_Coulomb_centimeter = constants.get("elementary charge") * constants.get(
        "Bohr radius") * constants.conversion_factor("m", "cm")

    numerator = 4.0 * np.pi**2 * N_A
    denominator = 3 * 1000 * np.log(10) * (4 * np.pi * e_0) * hbar * c

    return (numerator / denominator) * au_to_Coulomb_centimeter**2


def prefactor_ecd() -> float:
    r"""Prefactor for converting microscopic observable to decadic molar
    extinction coefficient in electronic circular dichroism.

    Notes
    -----
    This function implements the calculation of the following prefactor:

    .. math::

        k = \frac{16\pi^{2}N_{\mathrm{A}}}{3\times 1000\times \ln(10) (4 \pi \epsilon_{0}) n \hbar c^{2}}

    The prefactor is computed in SI units and then adjusted for the fact that
    we use atomic units to express microscopic observables: excitation energies
    and transition dipole moments.
    The refractive index :math:`n` is, in general, frequency-dependent. We
    assume it to be constant and equal to 1.

    """

    N_A = constants.get("Avogadro constant")
    c = constants.get("speed of light in vacuum")
    hbar = constants.get("Planck constant over 2 pi")
    e_0 = constants.get("electric constant")

    au_to_Coulomb_centimeter = constants.get("elementary charge") * constants.get(
        "Bohr radius") * constants.conversion_factor("m", "cm")
    au_to_Joule_inverse_Tesla = 2.0 * constants.get("Bohr magneton") * constants.conversion_factor("m", "cm")
    conversion = au_to_Coulomb_centimeter * au_to_Joule_inverse_Tesla

    numerator = 16.0 * np.pi**2 * N_A
    denominator = 3 * 1000 * np.log(10) * (4 * np.pi * e_0) * hbar * c**2

    return (numerator / denominator) * conversion


def spectrum(*,
             poles: Union[List[float], np.ndarray],
             residues: Union[List[float], np.ndarray],
             kind: str = "opa",
             lineshape: str = "gaussian",
             gamma: float = 0.2,
             npoints: int = 5000,
             out_units: str = "nm") -> Dict[str, np.ndarray]:
    r"""One-photon absorption (OPA) or electronic circular dichroism (ECD)
    spectra with phenomenological line broadening.

    This function gives arrays of values ready to be plotted as OPA spectrum:

    .. math::

       \varepsilon(\omega) =
          \frac{4\pi^{2}N_{\mathrm{A}}\omega}{3\times 1000\times \ln(10) (4 \pi \epsilon_{0}) n \hbar c}
          \sum_{i \rightarrow j}g_{ij}(\omega)|\mathbf{\mu}_{ij}|^{2}

    or ECD spectrum:

    .. math::

       \Delta\varepsilon(\omega) =
          \frac{16\pi^{2}N_{\mathrm{A}}\omega}{3\times 1000\times \ln(10) (4 \pi \epsilon_{0}) n \hbar c^{2}}
          \sum_{i \rightarrow j}g_{ij}(\omega)\Im(\mathbf{\mu}_{ij}\cdot\mathbf{m}_{ij})

    in macroscopic units of :math:`\mathrm{L}\cdot\mathrm{mol}^{-1}\cdot\mathrm{cm}^{-1}`.
    The lineshape function :math:`g_{ij}(\omega)` with phenomenological
    broadening :math:`\gamma` is used for the convolution of the infinitely
    narrow results from a linear response calculation.

    Parameters
    ----------
    poles
        Poles of the response function, i.e. the excitation energies.
        These are **expected** in atomic units of angular frequency.
    residues
        Residues of the linear response functions, i.e. transition dipole moments (OPA) and rotatory strengths (ECD).
        These are **expected** in atomic units.
    kind
        {"opa", "ecd"}
        Which kind of spectrum to generate, one-photon absorption ("opa") or electronic circular dichroism ("ecd").
        Default is `opa`.
    lineshape
        {"gaussian", "lorentzian"}
        The lineshape function to use in the fitting. Default is `gaussian`.
    gamma
        Full width at half maximum of the lineshape function.
        Default is 0.2 au of angular frequency.
        This value is **expected** in atomic units of angular frequency.
    npoints
        How many points to generate for the x axis. Default is 5000.
    out_units
        Units for the output array `x`, the x axis of the spectrum plot.
        Default is wavelengths in nanometers.
        Valid (and case-insensitive) values for the units are:

          - `au` atomic units of angular frequency
          - `Eh` atomic units of energy
          - `eV`
          - `nm`
          - `THz`

    Returns
    -------
    spectrum : Dict[str, numpy.ndarray]
        The fitted electronic absorption spectrum, with units for the x axis specified by the `out_units` parameter.
        This is a dictionary containing the convoluted (key: `convolution`) and the infinitely narrow spectra (key: `sticks`).

        .. code-block:: python

           {"convolution": {"x": np.ndarray, "y": np.ndarray},
            "sticks": {"poles": np.ndarray, "residues": np.ndarray}}

    Notes
    -----
    * Conversion of the broadening parameter :math:`\gamma`.
      The lineshape functions are formulated as functions of the angular frequency :math:`\omega`.
      When converting to other physical quantities, the broadening parameter has to be modified accordingly.
      If :math:`\gamma_{\omega}` is the chosen broadening parameter then:

        - Wavelength: :math:`gamma_{\lambda} = \frac{\lambda_{ij}^{2}}{2\pi c}\gamma_{\omega}`
        - Frequency: :math:`gamma_{\nu} = \frac{\gamma_{\omega}}{2\pi}`
        - Energy: :math:`gamma_{E} = \gamma_{\omega}\hbar`

    References
    ----------
    A. Rizzo, S. Coriani, K. Ruud, "Response Function Theory Computational Approaches to Linear and Nonlinear Optical Spectroscopy". In Computational Strategies for Spectroscopy.
    https://doi.org/10.1002/9781118008720.ch2
    """

    # Transmute inputs to np.ndarray
    if isinstance(poles, list):
        poles = np.array(poles)
    if isinstance(residues, list):
        residues = np.array(residues)
    # Validate input arrays
    if poles.shape != residues.shape:
        raise ValueError(f"Shapes of poles ({poles.shape}) and residues ({residues.shape}) vectors do not match!")

    # Validate kind of spectrum
    kind = kind.lower()
    valid_kinds = ["opa", "ecd"]
    if kind not in valid_kinds:
        raise ValueError(f"Spectrum kind {kind} not among recognized ({valid_kinds})")

    # Validate output units
    out_units = out_units.lower()
    valid_out_units = ["au", "eh", "ev", "nm", "thz"]
    if out_units not in valid_out_units:
        raise ValueError(f"Output units {out_units} not among recognized ({valid_out_units})")

    c = constants.get("speed of light in vacuum")
    c_nm = c * constants.conversion_factor("m", "nm")
    hbar = constants.get("Planck constant over 2 pi")
    h = constants.get("Planck constant")
    Eh = constants.get("Hartree energy")
    au_to_nm = 2.0 * np.pi * c_nm * hbar / Eh
    au_to_THz = (Eh / h) * constants.conversion_factor("Hz", "THz")
    au_to_eV = constants.get("Hartree energy in eV")

    converters = {
        "au": lambda x: x,  # Angular frequency in atomic units
        "eh": lambda x: x,  # Energy in atomic units
        "ev": lambda x: x * au_to_eV,  # Energy in electronvolts
        "nm": lambda x: au_to_nm / x,  # Wavelength in nanometers
        "thz": lambda x: x * au_to_THz,  # Frequency in terahertz
    }

    # Perform conversion of poles from au of angular frequency to output units
    poles = converters[out_units](poles)

    # Broadening functions
    gammas = {
        "au": lambda x_0: gamma,  # Angular frequency in atomic units
        "eh": lambda x_0: gamma,  # Energy in atomic units
        "ev": lambda x_0: gamma * au_to_eV,  # Energy in electronvolts
        "nm": lambda x_0: ((x_0**2 * gamma * (Eh / hbar)) / (2 * np.pi * c_nm)),  # Wavelength in nanometers
        "thz": lambda x_0: gamma * au_to_THz,  # Frequency in terahertz
    }

    # Generate x axis
    # Add a fifth of the range on each side
    expand_side = (np.max(poles) - np.min(poles)) / 5
    x = np.linspace(np.min(poles) - expand_side, np.max(poles) + expand_side, npoints)

    # Validate lineshape
    lineshape = lineshape.lower()
    valid_lineshapes = ["gaussian", "lorentzian"]
    if lineshape not in valid_lineshapes:
        raise ValueError(f"Lineshape {lineshape} not among recognized ({valid_lineshapes})")

    # Obtain lineshape function
    shape = Gaussian(x, gammas[out_units]) if lineshape == "gaussian" else Lorentzian(x, gammas[out_units])

    # Generate y axis, i.e. molar decadic absorption coefficient
    prefactor = prefactor_opa() if kind == "opa" else prefactor_ecd()
    transform_residue = (lambda x: x**2) if kind == "opa" else (lambda x: x)
    y = prefactor * x * np.sum([transform_residue(r) * shape.lineshape(p) for p, r in zip(poles, residues)], axis=0)

    # Generate sticks
    sticks = prefactor * np.array([p * transform_residue(r) * shape.maximum(p) for p, r in zip(poles, residues)])

    return {"convolution": {"x": x, "y": y}, "sticks": {"poles": poles, "residues": sticks}}
