# ⚡ ATPase-EP-Analyzer ⚡  
*A Nextflow workflow for ATP synthase electrostatic potential (ESP) analysis*

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A525.04.6-brightgreen?logo=nextflow)](https://www.nextflow.io/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue)](LICENSE)

---

## 🔬 Overview
**ATPase-EP-Analyzer** is a **Nextflow-based workflow** that automates the full pipeline of  
electrostatic potential (ESP) analysis for **F-type ATP synthase enzymes**:

1. **Alignment**: Custom rotor-focused structural alignment of PDBx/mmCIF input using pseudoatom *L keys* in PyMOL.  
2. **Fixing**: Fixing aligned PDB structures using PDBFixer.
3. **PQR Generation**: Assigning charges and protonation states to fixed PDB structures via PDB2PQR with configurable pH.  
4. **APBS Calculation**: Constructing electrostatic potential (ESP) maps with APBS.  
5. **Analysis**: Generating ESP “CT-scans,” angular averaging, and electric field estimation.  
6. **Visualization**: Compiling publication-ready figures & summary statistics.

---

## 📊 Workflow
<p align="center">
  <img src="https://raw.githubusercontent.com/Orbitaliant/ATPase-EP-Analyzer/main/assets/workflow_diagram.png" width="600" alt="Workflow diagram"/>
</p>

*(Tip: replace with a Nextflow DAG screenshot or schematic of your workflow.)*

---

## 🚀 Quickstart

### Requirements
- Java 17+
- [Nextflow ≥ 25.04.6.5954](https://www.nextflow.io/)
- Conda/Mamba

### Run Example
```bash
nextflow run ATPase-EP-Analyzer.nf --pH 7.4 --distmin 5.0 --distmax 15.0 --cif_input_dir ./cif_input
```

---

## 📥 Inputs
- **CIF/PDB files** in `./cif_input/` (or custom path via `--cif_input_dir`).  
- Parameters:
  - `--pH` (default `7.0`) → protonation state for PDB2PQR  
  - `--distmin` / `--distmax` → radial analysis window in Å  
  - `--openmm_platform` (default `CUDA`) → molecular mechanics platform  

---

## 📤 Outputs
- `results/atpase/<pdbid>/` → ATP synthase processed structure files, APBS maps, plots  
- `results/rotor/<pdbid>/`  → Rotor processed structure files, APBS maps, plots  
- `reports/` → Nextflow run reports, trace, timeline  
- `work/`   → Nextflow work directory (cache)  

---

## 🛠️ Installation
Clone the repository:
```bash
git clone https://github.com/Orbitaliant/ATPase-EP-Analyzer.git
cd ATPase-EP-Analyzer
```

Test the workflow with an example:
```bash
nextflow run ATPase-EP-Analyzer.nf --cif_input_dir assets/example_cif
```

---

## 📖 Citation
If you use **ATPase-EP-Analyzer**, please cite:

```yaml
@software{Matar2025ATPaseEP,
  author  = {Islam K. Matar},
  title   = {ATPase-EP-Analyzer: Nextflow workflow for ATP synthase ESP analysis},
  version = {v0.1.0},
  year    = {2025},
  url     = {https://github.com/Orbitaliant/ATPase-EP-Analyzer}
}
```

See [`CITATION.cff`](CITATION.cff) for more details.

---

## 🤝 Contributing
Contributions are welcome!  
Please:
- Open an [issue](https://github.com/Orbitaliant/ATPase-EP-Analyzer/issues) for bug reports & feature requests.  
- Use pull requests for improvements (documentation, code, tests).  

---

## 📜 License
Licensed under the **Apache-2.0 License**.  
See the [LICENSE](LICENSE) file for details.

---

<p align="center">
✨ Developed with Nextflow for scalable & reproducible ATP synthase research ✨
</p>
