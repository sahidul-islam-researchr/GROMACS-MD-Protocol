# Alzheimer's Drug Discovery: Multi-Target Phytochemical Inhibitors

## Overview
This repository contains the computational workflow and methodology for the identification of potent phytochemical inhibitors against Acetylcholinesterase (AChE) and Butyrylcholinesterase (BChE) for Alzheimer's Disease (AD) therapy. The study utilizes an integrated *in silico* approach combining molecular docking, ADMET profiling, and 100 ns molecular dynamics simulations.

## Research Abstract
Alzheimer's disease is addressed by targeting cholinergic receptors AChE (PDB: 7AIY) and BChE (PDB: 7E3H). This project screens 18 phytoconstituents from 11 medicinal plant families to identify multi-target inhibitors with superior pharmacological profiles compared to existing FDA-approved treatments.
![RMSD Analysis of Protein-Ligand Complex](images/filename.png)

## Computational Workflow & Tools
The integrated pipeline covers the following stages:

- **Target Preparation:** PDB structure processing for 7AIY and 7E3H.
- **Molecular Docking:** Performed using **AutoDock Vina** to evaluate binding affinities.
- **ADMET & Toxicity Profiling:** Analyzed using **SwissADME** and **ProTox 3.0**.
- **Molecular Dynamics (MD) Simulations:** 100 ns all-atom simulations conducted via **GROMACS**.
- **Thermodynamic Analysis:** Binding free energy calculations performed using **MM-GBSA**.

## Key Features
- **Multi-Target approach:** Dual inhibition strategy for AChE and BChE.
- **Dynamic Stability:** Analysis of protein-ligand stability over 100ns trajectory.
- **Pharmacological Insight:** Evaluation of therapeutic efficiency, stability, and safety of natural compounds.

## Citation
If you find this work helpful for your research, please refer to the documentation or contact for further inquiries regarding the methodology.

*Maintained by: Sahidul Islam *
