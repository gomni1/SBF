# SBF Verification Suite (Python)

This directory contains the source code used to generate the falsifiable predictions of the Single Bulk Framework.

## 🛠️ Repository Files
* **`sbf_verification_suite.py`**: The Master Execution Script. Runs the full simulation and generates compliance reports.
* **`sbf_topology.py`**: The Topological Kernel (Phase 2). Simulations of Dark Matter/Void scaling.
* **`sbf_core.py`**: The Physics Engine. Handles Lepton Mass Hierarchy and Neutron calculations.
* **`Topological_Nyquist_Limit.py`**
* **`sbf_engine_v3.py`**

* "The N-body solution is verified through the Topological Nyquist Limit script, which defines the geometric cutoff for lattice interactions. The mechanical execution of this limit is handled by the calculate_torsion_displacement function in sbf_engine_v3.py, and the stability of the resulting scalar sum is validated in the Sector 3 checks of the sbf_verification_suite.py."

## 🚀 How to Run the Code

### Prerequisites
* Python 3.8+
* NumPy
* Matplotlib

### Installation & Execution
```bash
# Clone the repository
git clone [https://github.com/gomni1/SBF.git](https://github.com/gomni1/SBF.git)

# Navigate to this directory
cd "SBF/SBF PYTHON"

# Install dependencies
pip install -r requirements.txt

# Run the Verification Suite
python sbf_verification_suite.py
