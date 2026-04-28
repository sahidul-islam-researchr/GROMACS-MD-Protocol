# GROMACS MD Simulation & Analysis Protocol for Protein-Ligand System

This repository contains all commands and workflows to run a molecular dynamics (MD) simulation and subsequent analysis using **GROMACS (CHARMM36 + TIP3P)**.  
The protocol is designed for publication-ready results.

## 📁 Repository Structure
- `README.md` – Main workflow (copy/paste commands step by step)
- `results/` – Analysis outputs (`.xvg`, `.png`, `.csv`, etc.) – you will fill this after running

## 🧪 Prerequisites
- GROMACS (2020+) sourced: `source /usr/local/gromacs/bin/GMXRC`
- Conda environments: `pdbfixer`, `gmxMMPBSA`, `md_analysis`
- SwissParam web server for ligand topology
- Grace (`xmgrace`) for plotting

---

## 1️⃣ Simulation Workflow

### 1.1 Prepare Protein and Ligand Files

```bash
# Fix protein with pdbfixer
conda activate pdbfixer
pdbfixer REC.pdb --output fixed_REC.pdb

# Generate protein topology
gmx pdb2gmx -f fixed_REC.pdb -ignh -ff charmm36 -water tip3p

# For ligand: upload LIG.mol2 to SwissParam → download .zip → extract LIG.itp, LIG.gro, LIG.pdb
