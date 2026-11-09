# Load derivation script for a hot-staging metallic truss interstage

This repository contains the script that computes every number in Section 3 (Load Derivations) of the paper

> C. Hartsink, \*Failure-Mode-Driven FEA-Based Design Optimization of Hot-Staging Metallic Lattice Truss Interstage Structures Under the Reusable-Booster Flight Cycle\*, International Journal of Structural Stability and Dynamics (submitted).

The paper studies the failure-mode-driven design of a metallic lattice-truss hot-staging interstage, with Starship V3 as the reference vehicle. The script derives the mechanical loads at maximum dynamic pressure (max-q) and just before main engine cutoff (pre-MECO), and the convective heat load of the upper-stage engine plume during hot-staging.

## Contents

|Path|Description|
|-|-|
|`interstage\_loads.py`|The load derivation script. It runs top to bottom and prints its results in the order of the subsections of Section 3.|
|`output/interstage\_loads\_output.txt`|The unedited console output of the run reported in the paper.|
|`cea\_inputs/`|The three input files for NASA's Chemical Equilibrium with Applications (CEA) program listed in Appendix A of the paper: the rocket problem for the nozzle exit states (A.1), the assigned enthalpy and pressure problem for the stagnation state (A.2), and the assigned temperature and pressure problem for the wall and reference states (A.3).|
|`requirements.txt`|The Python package versions of the reported run.|

## Running the script

The reported run used Python 3.13 with the package versions in `requirements.txt`:

```
python -m venv .venv
source .venv/bin/activate        # on Windows: .venv\\Scripts\\activate
pip install -r requirements.txt
python interstage\_loads.py
```

Its output should match `output/interstage\_loads\_output.txt`. RocketCEA may print additional lines about its installation directory when it is first imported; these are not part of the results.

* **RocketCEA** calls the CEA Fortran code (CEA2) for the combustion and nozzle expansion of the Raptor engines.
* **ambiance** provides the 1976 U.S. Standard Atmosphere for the ascent heating estimate of Section 3.4.
* **Cantera** is optional. It is only used to cross-check the stagnation temperature against an independent equilibrium calculation (Section 3.7); without it, the script skips that check.

The values of the CEA problem of Appendix A.3, which RocketCEA does not solve, were computed with the NASA CEARUN web interface and are entered directly in the script. The input files in `cea\_inputs/` reproduce those runs in CEARUN or the stand-alone CEA program.

## Citation

If you use this script, please cite the paper above. Citation metadata for the script itself is in `CITATION.cff`.

## License

The script is released under the MIT License (see `LICENSE`). It imports RocketCEA, which is licensed under the GNU General Public License v3; RocketCEA is not included in this repository and is installed separately.

