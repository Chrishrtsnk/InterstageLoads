# Load derivation script for a hot-staging metallic truss interstage

This repository contains the script that computes every number in Section 3 (Load Derivations) and Appendix A.2 of the paper

> C. Hartsink and S. Verma, *Failure-Mode-Driven FEA-Based Design Optimization of Hot-Staging Metallic Lattice Truss Interstage Structures Under the Reusable-Booster Flight Cycle*, International Journal of Structural Stability and Dynamics (submitted).

The paper studies the failure-mode-driven design of a metallic lattice-truss hot-staging interstage, with Starship V3 as the reference vehicle. The script derives the mechanical loads at maximum dynamic pressure (max-q) and just before main engine cutoff (pre-MECO), the ascent heating before separation, and the convective heat load of the upper-stage engine plume during hot-staging, and it repeats the validation and sensitivity checks of Appendix A.2.6, including the check of the pre-MECO load factor against flight telemetry.

## Contents

| Path | Description |
|---|---|
| `interstage_loads.py` | The load derivation script. It runs top to bottom; each block of output is labeled with the section of the paper and the appendix that contains its derivation, e.g. `[3.2 \| A.2.2]`. |
| `output/interstage_loads_output.txt` | The unedited console output of the run reported in the paper. |
| `cea_inputs/` | The three input files for NASA's Chemical Equilibrium with Applications (CEA) program listed in Appendix A.2.7 of the paper: the rocket problem for the nozzle exit states (`A1_rocket_exit_states.inp`), the assigned enthalpy and pressure problem for the stagnation state (`A2_hp_stagnation_state.inp`) and the assigned temperature and pressure problem for the wall and reference states (`A3_tp_wall_and_reference_states.inp`). |
| `requirements.txt` | The Python package versions of the reported run. |

## Where each result appears in the paper

| Output label | Paper |
|---|---|
| `[3.2 \| A.2.2]` | Axial compression: Section 3.2, Appendix A.2.2 |
| `[3.2 \| A.2.3]` | Shear force and bending moment: Section 3.2, Appendix A.2.3 |
| `[3.3 \| A.2.4]` | Pre-staging thermal state and ascent heating: Section 3.3, Appendix A.2.4 |
| `[3.3 \| A.2.5]` | Plume impingement chain and equivalent heat transfer coefficient: Section 3.3, Appendix A.2.5 |
| `[3.4]` | Limit, yield and ultimate loads: Section 3.4, Table 3.1 |
| `[A.2.6]` | Validation, sensitivity (Table A.2.3) and telemetry check: Appendix A.2.6 |

## Running the script

The reported run used Python 3.13 with the package versions in `requirements.txt`:

```
python -m venv .venv
source .venv/bin/activate        # on Windows: .venv\Scripts\activate
pip install -r requirements.txt
python interstage_loads.py
```

Its output should match `output/interstage_loads_output.txt`. RocketCEA may print additional lines about its installation directory when it is first imported; these are not part of the results. Installing RocketCEA from source requires a Fortran compiler (for example gfortran).

* **RocketCEA** calls the CEA Fortran code (CEA2) for the combustion and nozzle expansion of the Raptor engines.
* **ambiance** provides the 1976 U.S. Standard Atmosphere for the ascent heating estimate of Appendix A.2.4. The max-q loads use the tabulated values at 10 km geopotential altitude.
* **Cantera** is optional. It is only used to cross-check the stagnation temperature against an independent equilibrium calculation (Appendix A.2.6); without it, the script skips that check.

The values of the CEA problem `A3_tp_wall_and_reference_states.inp`, which RocketCEA does not solve, were computed with the NASA CEARUN web interface and are entered directly in the script. The input files in `cea_inputs/` reproduce those runs in CEARUN or the stand-alone CEA program.

## Citation

If you use this script, please cite the paper above. Citation metadata for the script itself is in `CITATION.cff`.

## License

The script is released under the MIT License (see `LICENSE`). It imports RocketCEA, which is licensed under the GNU General Public License v3; RocketCEA is not included in this repository and is installed separately.
