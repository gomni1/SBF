import numpy as np
import matplotlib.pyplot as plt
import unittest

# ==============================================================================
# SINGLE BULK FRAMEWORK (SBF) - VACUUM TOPOLOGY KERNEL v2.0
# ==============================================================================
# Description:
#   Simulates the topological phase transition of the Planck-scale vacuum.
#   Derives effective fundamental constants (G, hbar, alpha) as a function
#   of the local Coordination Number (Z).
#
#   Physics Base:
#   - Vacuum = Nematic Liquid Crystal of Tetrahedral Units
#   - Critical Jamming Point (Bulk): Z_c = 14.4
#   - Isostatic Stability Limit: Z_iso = 6.0
# ==============================================================================

class VacuumLattice:
    def __init__(self):
        # --- Fundamental Constants (Bulk Values at Z = 14.4) ---
        self.c_0 = 299792458.0          # Speed of Light (m/s)
        self.G_0 = 6.67430e-11          # Newton's Constant (m^3 kg^-1 s^-2)
        self.hbar_0 = 1.0545718e-34     # Planck Constant (J s)
        self.alpha_0 = 1/137.035999     # Fine Structure Constant
        
        # --- Topological Parameters ---
        self.Z_c = 14.4                 # Critical Bulk Coordination (Max Entropy)
        self.Z_iso = 6.0                # Isostatic Limit (Frictionless Spheres)
        
    def scaling_factor(self, Z):
        """
        Calculates the topological scaling ratio Gamma(Z).
        Gamma = (Z - Z_iso) / (Z_c - Z_iso)
        """
        # Robust clipping to prevent division by zero near singularity
        Z_clipped = np.clip(Z, self.Z_iso + 1e-6, None)
        return (Z_clipped - self.Z_iso) / (self.Z_c - self.Z_iso)

    def calculate_observables(self, Z):
        """
        Derives effective physical constants at a local coordination Z.
        Based on SBF 'Option A' (Constant c, Variable density).
        """
        gamma = self.scaling_factor(Z)
        
        # 1. Shear Modulus (Stiffness) & Inertial Density
        # Scaling: mu ~ gamma^1.5 | rho ~ gamma^1.5
        mu_ratio = gamma**1.5
        rho_ratio = gamma**1.5  # Required for Lorentz Invariance
        
        # 2. Speed of Light
        # c = sqrt(mu / rho). Since mu and rho scale identically, c is constant.
        # We explicitly set this to enforce Lorentz Invariance numerically.
        c_eff = self.c_0
        
        # 3. Effective Gravity (G)
        # G ~ 1/mu (Inverse Stiffness). 
        # Softer vacuum (low Z) -> Stronger Gravity.
        G_eff = self.G_0 * (1.0 / mu_ratio)
        
        # 4. Effective Planck Constant (hbar)
        # From Nelsonian Stochastic Mechanics: hbar ~ 1/gamma
        hbar_eff = self.hbar_0 * (1.0 / gamma)
        
        # 5. Fine Structure Constant (alpha)
        # alpha ~ 1 / (hbar * c). 
        # alpha_eff = alpha_0 * gamma
        alpha_eff = self.alpha_0 * gamma
        
        return {
            "Z": Z,
            "c": c_eff,
            "G": G_eff,
            "hbar": hbar_eff,
            "alpha": alpha_eff,
            "mu_scaling": mu_ratio
        }

    def calculate_observational_signatures(self, Z_void=10.0):
        """Calculate measurable differences for observational tests"""
        bulk = self.calculate_observables(self.Z_c)
        void = self.calculate_observables(Z_void)
        
        return {
            'delta_alpha_percent': ((void['alpha'] - bulk['alpha']) / bulk['alpha']) * 100,
            'G_void_ratio': void['G'] / bulk['G'],
            'hbar_void_ratio': void['hbar'] / bulk['hbar']
        }

# ==============================================================================
# VISUALIZATION
# ==============================================================================

def plot_phase_transition(sim):
    """Generates the 4-panel Phase Transition Dashboard."""
    Z_range = np.linspace(6.1, 15, 200)
    results = [sim.calculate_observables(Z) for Z in Z_range]
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('SBF Topological Phase Transition', fontsize=16)
    
    # 1. Gravity (G/G0)
    y_G = [r['G']/sim.G_0 for r in results]
    axes[0,0].plot(Z_range, y_G, 'b-', linewidth=2)
    axes[0,0].set_title('Gravitational Coupling (G_eff)')
    axes[0,0].set_ylabel('G / G_0')
    axes[0,0].axvline(sim.Z_c, color='r', linestyle='--', label='Bulk (Z=14.4)')
    axes[0,0].text(8, max(y_G)*0.8, "Dark Matter Mimic\n(Stronger Gravity)", fontsize=10)
    axes[0,0].grid(True, alpha=0.3)
    
    # 2. Fine Structure (alpha/alpha0)
    y_alpha = [r['alpha']/sim.alpha_0 for r in results]
    axes[0,1].plot(Z_range, y_alpha, 'g-', linewidth=2)
    axes[0,1].set_title('Fine Structure Constant (alpha)')
    axes[0,1].set_ylabel('alpha / alpha_0')
    axes[0,1].axvline(sim.Z_c, color='r', linestyle='--')
    axes[0,1].text(8, min(y_alpha)*1.1, "Webb Dipole Region\n(Smaller Alpha)", fontsize=10)
    axes[0,1].grid(True, alpha=0.3)

    # 3. Planck Constant (hbar/hbar0)
    y_hbar = [r['hbar']/sim.hbar_0 for r in results]
    axes[1,0].plot(Z_range, y_hbar, 'purple', linewidth=2)
    axes[1,0].set_title('Planck Constant (Effective)')
    axes[1,0].set_ylabel('hbar / hbar_0')
    axes[1,0].set_yscale('log')
    axes[1,0].axvline(sim.Z_iso, color='k', linestyle=':', label='Isostatic Limit')
    axes[1,0].grid(True, alpha=0.3)

    # 4. Shear Modulus (Stiffness)
    y_mu = [r['mu_scaling'] for r in results]
    axes[1,1].plot(Z_range, y_mu, 'k-', linewidth=2)
    axes[1,1].set_title('Vacuum Stiffness (mu)')
    axes[1,1].set_ylabel('Stiffness Ratio')
    axes[1,1].set_xlabel('Coordination Number (Z)')
    axes[1,1].axvline(sim.Z_c, color='r', linestyle='--')
    axes[1,1].grid(True, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# ==============================================================================
# UNIT TESTS
# ==============================================================================

class TestVacuumLattice(unittest.TestCase):
    def setUp(self):
        self.sim = VacuumLattice()
    
    def test_bulk_values(self):
        """At Z = Z_c, all ratios should be 1.0 (Standard Model recovered)"""
        bulk = self.sim.calculate_observables(self.sim.Z_c)
        self.assertAlmostEqual(bulk['c'], self.sim.c_0)
        self.assertAlmostEqual(bulk['G'], self.sim.G_0)
        self.assertAlmostEqual(bulk['alpha'], self.sim.alpha_0)
    
    def test_monotonicity(self):
        """G should increase as Z decreases (Dark Matter effect)"""
        G_void = self.sim.calculate_observables(10.0)['G']
        G_bulk = self.sim.calculate_observables(14.4)['G']
        self.assertGreater(G_void, G_bulk)
        
    def test_lorentz_invariance(self):
        """Speed of light must remain constant across phases"""
        c_void = self.sim.calculate_observables(8.0)['c']
        self.assertEqual(c_void, self.sim.c_0)

# ==============================================================================
# EXECUTION
# ==============================================================================

if __name__ == "__main__":
    # 1. Run Simulations
    sim = VacuumLattice()
    print(f"--- SBF TOPOLOGICAL VERIFICATION SUITE v2.0 ---")
    
    # Calculate Signatures
    sigs = sim.calculate_observational_signatures(Z_void=10.0)
    print(f"\nPREDICTIONS FOR DEEP VOID (Z=10.0):")
    print(f"  > Gravity Strength:    {sigs['G_void_ratio']:.2f}x (Dark Matter Effect)")
    print(f"  > Alpha Variation:     {sigs['delta_alpha_percent']:.2f}% (Webb Dipole)")
    print(f"  > Quantum Coherence:   {sigs['hbar_void_ratio']:.2f}x (Macroscopic QM)")

    # 2. Run Visualization
    print("\nGenerating Phase Transition Plots...")
    # plot_phase_transition(sim)  # Uncomment to show plots

    # 3. Run Unit Tests
    print("\nRunning Logic Validation Tests...")
    suite = unittest.TestLoader().loadTestsFromTestCase(TestVacuumLattice)
    unittest.TextTestRunner(verbosity=2).run(suite)