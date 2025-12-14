"""
SINGLE BULK FRAMEWORK (SBF) - CORE PHYSICS ENGINE - v5.0 (Final Release)
Author: Glenn Millar
Date: December 2025

Description:
    The mathematical kernel of the Single Bulk Framework. Defines geometric axioms,
    physical constants, and constitutive laws. Includes inverse calculators and
    experimental comparison tools.

License:
    The MIT License (MIT)
    Copyright (c) 2025 Glenn Millar

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
"""

import math
from dataclasses import dataclass, field
from typing import ClassVar, Dict, Tuple, Optional

class SBFError(Exception):
    """Custom exception for SBF validation failures."""
    pass

class SBFConstants:
    """Class holding standard, non-configurable physical constants (SI Units)."""
    C: ClassVar[float] = 2.99792458e8        # Speed of light [m/s]
    HBAR: ClassVar[float] = 1.0545718e-34    # Reduced Planck constant [J*s]
    MU_0: ClassVar[float] = 1.25663706e-6    # Vacuum permeability [N/A^2]
    G_N: ClassVar[float] = 6.67430e-11       # Gravitational constant [m^3 kg^-1 s^-2]
    
    # Derived Planck Length (L_P = sqrt(hbar * G / c^3)) [Source: 1703]
    L_P: ClassVar[float] = math.sqrt((HBAR * G_N) / (C**3)) 
    
    # Cosmological Rigidity Floor (Dark Energy Density approximation) [Source: 1288]
    RHO_LAMBDA_FLOOR: ClassVar[float] = 1.0e-9 # Pa (~1e-9 J/m^3)
    
    # Constants for QCD Renormalization Group Flow [Source: 1045]
    ALPHA_PLANCK: ClassVar[float] = 0.0215 
    
    # Conversion factors
    CONV_GEV2: ClassVar[float] = 3.89379e-32    # (hbar*c)^2 in [GeV^2 * m^2]
    HBAR_C_GEV_FM: ClassVar[float] = 0.197327   # hbar*c in GeV*fm

    @classmethod
    def calculate_planck_pressure(cls) -> float:
        """Calculates Planck pressure P_P = c^7/(hbar*G^2)."""
        c7 = cls.C**7
        return c7 / (cls.HBAR * cls.G_N**2)


@dataclass
class SBFConfiguration:
    """
    User-configurable inputs derived from fundamental axioms and cosmology.
    """
    # Axiomatic Input: Coordination Number (Z) [Source: 1659]
    Z_BERNAL: float = 14.39 
    
    # Phenomenological Input: Electron Mass Baseline [Source: 1658]
    M_ELECTRON_MEV: float = 0.510998 
    
    # Cosmological Input: Hubble Radius [Source: 1682]
    R_HUBBLE: float = 1.4e26 

    # Empirical Input: Dilatancy Coefficient (beta) [Source: 1269, 1303]
    BETA_DILATANCY: float = 1.0e-4

    def __post_init__(self):
        """Input validation upon initialization."""
        if not (12.0 < self.Z_BERNAL < 16.0):
            raise SBFError(f"Z_BERNAL ({self.Z_BERNAL}) must be within the physically relevant jamming range (12 to 16).")
        if self.M_ELECTRON_MEV <= 0:
            raise SBFError("Electron mass must be positive.")
        if self.R_HUBBLE <= 1.0e20:
            raise SBFError("Hubble radius must be a large, macroscopic length scale.")
        if not (1.0e-5 < self.BETA_DILATANCY < 1.0):
            raise SBFError("Dilatancy coefficient is expected to be a small positive number.")
            
        # Consistency Check for Hubble Constant (H0 ~ c / R_H)
        h0_si = SBFConstants.C / self.R_HUBBLE
        # Convert to km/s/Mpc (1 Mpc = 3.086e22 m)
        h0_kms_mpc = h0_si * (3.086e22 / 1000.0)
        # We don't raise an error here, but it's a useful check for the user
        # Expected ~70 km/s/Mpc


class SBFModel:
    """
    The Single Bulk Framework Physics Engine.
    Requires an SBFConfiguration object to run.
    """
    
    def __init__(self, config: SBFConfiguration):
        self.config = config
        self.Z = self.config.Z_BERNAL
        self.M_e = self.config.M_ELECTRON_MEV
        self.xi = 1.0 / (2.0 * math.pi) # Curvature penalty [Source: 1662]
        self.constants = SBFConstants()

    # =========================================================================
    # SECTOR 1: MATTER (LEPTONS)
    # =========================================================================
    
    def calculate_lepton_mass(self, N: int) -> float:
        """
        [Section 3.2] Derives Lepton mass (MeV) via Topological Knot Scaling.
        Formula: M_N = M_e * Z^(N-3) * [1 + xi * correction]
        """
        if N not in {3, 5, 6, 7}:
            raise SBFError(f"Knot complexity N={N} is not a valid Lepton knot (expected 3, 5, 6, or 7).")
        
        if N == 3:
            return self.M_e
        
        scale_factor = self.Z**(N - 3)
        
        # Steric Crowding / Curvature Penalty [Section 3.2.3]
        stiffness_correction = 1.0
        if N >= 6:
            # Assumes linear scaling of penalty multiplier as per manuscript (1*xi for N=6, 2*xi for N=7)
            penalty_multiplier = float(N - 5) 
            stiffness_correction = 1.0 + (penalty_multiplier * self.xi)
            
        return self.M_e * scale_factor * stiffness_correction

    def compare_lepton_masses(self) -> Dict[int, Tuple[float, Optional[float]]]:
        """Return predicted vs experimental lepton masses [Source: 1695]."""
        experimental = {
            3: 0.510998,  # Electron
            5: 105.658,   # Muon
            6: 1776.86,   # Tau
            7: None       # 4th gen (prediction)
        }
        return {N: (self.calculate_lepton_mass(N), experimental.get(N)) 
                for N in [3, 5, 6, 7]}

    # =========================================================================
    # SECTOR 2: ELECTROMAGNETISM
    # =========================================================================

    def calculate_fine_structure_inverse(self) -> float:
        """
        [Section 7.4.3] Derives inverse Fine Structure Constant (alpha^-1).
        Formula: alpha^-1 = (2/3) * Z^2
        """
        f_screen = 2.0 / 3.0
        return f_screen * (self.Z ** 2)

    # =========================================================================
    # SECTOR 3: COSMOLOGY (DARK ENERGY)
    # =========================================================================

    def calculate_dark_energy_density(self) -> float:
        """
        [Section 4.4] Derives Vacuum Energy Density (J/m^3) via Holographic Scaling.
        Formula: rho_lambda = (hbar * c) / (L_P^2 * R_H^2)
        """
        numerator = self.constants.HBAR * self.constants.C
        denominator = (self.constants.L_P**2) * (self.config.R_HUBBLE**2)
        return numerator / denominator

    # =========================================================================
    # SECTOR 4: ANOMALIES (NEUTRON LIFETIME)
    # =========================================================================

    def calculate_neutron_lifetime_shift(self, B_field_tesla: float = 1.0) -> float:
        """
        [Section 11.2.6.1] Calculates fractional shift in neutron decay rate.
        Formula: dGamma/Gamma = beta * ln(sigma_ext / rho_lambda)
        """
        if B_field_tesla <= 0:
            raise SBFError("Magnetic field strength must be positive.")
            
        sigma_ext = (B_field_tesla**2) / (2 * self.constants.MU_0)
        rho_lambda_floor = self.constants.RHO_LAMBDA_FLOOR 
        
        ratio = sigma_ext / rho_lambda_floor
        if ratio <= 1.0:
            return 0.0
            
        log_term = math.log(ratio)
        return self.config.BETA_DILATANCY * log_term

    def calculate_dilatancy_from_observation(self, observed_shift: float = 0.009, 
                                             B_field_tesla: float = 1.0) -> float:
        """
        Infer required beta (dilatancy) from observed neutron lifetime shift.
        Useful for checking the 'running' of beta between gravity and magnetic regimes.
        """
        sigma_ext = (B_field_tesla**2) / (2 * self.constants.MU_0)
        ratio = sigma_ext / self.constants.RHO_LAMBDA_FLOOR
        if ratio <= 1.0: return 0.0
        return observed_shift / math.log(ratio)

    # =========================================================================
    # SECTOR 5: STRONG FORCE (QCD)
    # =========================================================================

    def calculate_qcd_string_tension(self, L_target: float = 1.0e-15) -> float:
        """
        [Section 9] Derives QCD string tension (GeV^2).
        Formula: alpha(L) = alpha_P / (1 - alpha_P * ln(L/L_P))
        """
        if L_target <= self.constants.L_P:
            raise SBFError("Target length must be > Planck length.")
            
        alpha_P = self.constants.ALPHA_PLANCK
        ratio = L_target / self.constants.L_P
        log_term = math.log(ratio)
        
        # Denominator of RG flow equation
        denominator = 1.0 - (alpha_P * log_term)
        
        # Robust Pole Detection (DeepSeek recommendation)
        # Using relative tolerance check for stability near singularity
        if abs(denominator / (alpha_P * log_term)) < 1e-12:
            return float('inf') 
            
        alpha_L = alpha_P / denominator
        
        # Convert to Tension kappa in GeV^2
        kappa_geometric = alpha_L / (L_target**2) # m^-2
        return kappa_geometric * self.constants.CONV_GEV2

    def calculate_qcd_string_tension_gev_per_fm(self, L_target: float = 1.0e-15) -> float:
        """
        Returns string tension in the more common unit [GeV/fm].
        """
        kappa_GeV2 = self.calculate_qcd_string_tension(L_target)
        if kappa_GeV2 == float('inf'):
            return float('inf')
        # Convert GeV^2 -> GeV/fm
        return math.sqrt(kappa_GeV2) / self.constants.HBAR_C_GEV_FM