import numpy as np
import matplotlib.pyplot as plt

# --- SBF CONSTANTS ---
M_glueball = 1.89  # GeV (Scalar Glueball Mass / Yield Point)
alpha_inv = 137.036 # Network Stiffness
c = 1.0 # Natural units

# --- DERIVED LIMITS ---
# Static Mass Ceiling: M_max = M_0 * (alpha^-1)^2
M_max_TeV = (M_glueball * alpha_inv**2) / 1000.0

# Dynamic Bandwidth Limit (The Knee): E_knee = M_max * alpha^-1
E_knee_PeV = (M_max_TeV * alpha_inv) / 1000.0

print(f"SBF PREDICTIONS:")
print(f"Mass Ceiling (M_max): {M_max_TeV:.2f} TeV")
print(f"Cosmic Ray Knee (E_knee): {E_knee_PeV:.2f} PeV")

# --- DISPERSION RELATION FUNCTION ---
# E^2 = p^2 + m^2 - eta * (E^4 / E_knee^2)
def sbf_dispersion(p_array, mass_amu, eta=1.0):
    # E_knee scales with nucleon count (A) for composite nuclei
    # A Proton has A=1, Iron has A=56
    scaled_knee = E_knee_PeV * mass_amu * 1e6 # Convert PeV to GeV
    
    # We solve E for each p (approximated for high energy E ~ p)
    # The penalty term is -eta * (p^4 / E_knee^2)
    E_sbf = []
    for p in p_array:
        # Standard Model Energy
        E_sm = np.sqrt(p**2 + (mass_amu * 0.938)**2) # Mass in GeV
        
        # SBF Correction factor (Topological Cherenkov Braking)
        # We apply a soft suppression as p approaches the limit
        suppression = 1.0 - eta * (p / scaled_knee)**3 
        
        if suppression < 0: suppression = 0 # Energy saturation
        
        E_sbf.append(E_sm * suppression)
    return np.array(E_sbf)

# --- PLOTTING DATA ---
# Momentum range: 100 TeV to 1000 PeV (Log scale)
p_range = np.logspace(5, 9, 500) # In GeV

# 1. Standard Model (Pure Linear)
E_sm = p_range 

# 2. SBF Proton (A=1)
E_proton = sbf_dispersion(p_range, mass_amu=1)

# 3. SBF Iron Nucleus (A=56)
E_iron = sbf_dispersion(p_range, mass_amu=56)

# --- THE PLOT ---
plt.figure(figsize=(10, 6))

# Plot Standard Model Reference
plt.loglog(p_range/1e6, E_sm/1e6, 'k--', label='Standard Model (Lorentz Invariant)', alpha=0.5)

# Plot SBF Proton
plt.loglog(p_range/1e6, E_proton/1e6, 'b-', linewidth=2, label=f'SBF Proton (Knee @ {E_knee_PeV:.2f} PeV)')

# Plot SBF Iron
plt.loglog(p_range/1e6, E_iron/1e6, 'r-', linewidth=2, label=f'SBF Iron (Knee @ {E_knee_PeV*56:.1f} PeV)')

# Formatting
plt.axvline(x=E_knee_PeV, color='b', linestyle=':', alpha=0.6)
plt.axvline(x=E_knee_PeV*56, color='r', linestyle=':', alpha=0.6)

plt.xlabel('Momentum (PeV)')
plt.ylabel('Energy (PeV)')
plt.title('SBF Predicted Deviation from Lorentz Invariance (The Knee)')
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.legend()
plt.text(0.2, 2.0, "SBF Vacuum Aliasing Region", rotation=45, fontsize=12, color='gray')

plt.savefig("SBF_Dispersion_Knee.png", dpi=300)
plt.show()