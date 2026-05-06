# Bacterial Mutation Bias Analysis

A bioinformatics pipeline for characterising substitution mutation bias across 14 bacterial species, developed as a final year project at the University of Bath (2025–2026).

The pipeline identifies singleton SNPs from whole-genome sequencing data, filters likely recombination events, builds substitution spectra normalised by genome composition, and classifies mutations by functional consequence using snpEff. Results are used to test whether mutation bias varies systematically across species with different ecological niches and GC contents.

---

## Species analysed

*Bordetella pertussis*, *Campylobacter jejuni*, *Escherichia coli*, *Haemophilus influenzae*, *Klebsiella pneumoniae*, *Listeria monocytogenes*, *Mycobacterium tuberculosis*, *Neisseria meningitidis*, *Pseudomonas aeruginosa*, *Salmonella typhimurium*, *Staphylococcus aureus*, *Staphylococcus epidermidis*, *Streptococcus agalactiae*, *Streptococcus pneumoniae*

---

## Repository contents

| File | Description |
|------|-------------|
| `Multi_species_analysis_14sp.ipynb` | Main analysis pipeline (Python). Steps 1–4: VCF filtering, singleton and recombination removal, mutation matrix construction, pyrimidine collapse, and genome composition normalisation. Steps 5–6: snpEff-based functional classification, trinucleotide context analysis, and CSV export for statistical modelling and figure generation. |
| `Figures.ipynb` | Figure generation (Python/matplotlib). Produces all figures from the CSV outputs of the analysis notebook. |
| `run_snpeff.sh` | Bash script for snpEff 5.4c installation, database setup, and VCF annotation for all 14 species. Must be run in the terminal between Steps 4 and 5 of the analysis notebook. |
| `poisson_regression_model.R` | Poisson regression model fitted in R. Tests for significant differences in mutation rate by substitution class, species, and coding/intergenic status. Exports results to CSV. |
| `model_results_14sp.csv` | Exported model coefficients, standard errors, z-values, raw p-values, and BH-adjusted p-values for all model terms. |
| `mutation_data_14sp.csv` | Per-species, per-substitution-class, per-category mutation counts and opportunity denominators. Input for the R model and figures notebook. |
| `mutation_data_category_opp.csv` | Per-species, per-category total mutation counts and opportunity denominators. Input for the dN/dS figure. |
| `trinucleotide_data_14sp.csv` | Per-species trinucleotide context counts, opportunities, and proportions across all 96 COSMIC contexts. |
| `gc_content_14sp.csv` | Genomic GC content per species. Input for the GC content vs C→T scatter figure. |

---

## Dependencies

**Python (analysis and figures notebooks)**
- Python 3.x
- pandas
- numpy
- matplotlib
- seaborn
- scipy

**Bash (snpEff pipeline)**
- conda
- snpEff 5.4c (installed via conda, see `run_snpeff.sh`)

---
## How to run

### 1. Run the analysis notebook (Steps 1–4)
Open `Multi_species_analysis_14sp.ipynb` in JupyterLab and run cells up to and including Step 4 (Kernel → Restart & Run All up to the normalisation step). This will:
- Copy raw Parsnp VCFs from the shared team directory into `vcfs/raw/`
- Apply singleton and recombination filters
- Build substitution matrices and collapse to 6 pyrimidine-convention classes
- Normalise by reference genome nucleotide composition

### 2. Run the snpEff annotation pipeline
Before running Step 5 of the analysis notebook, annotated VCFs must be generated for all 14 species. From the terminal on CLIMB:
```bash
cd /home/jovyan/shared-team/2025-masters-project/people/alison
bash run_snpeff.sh
```
This installs snpEff 5.4c, builds custom databases from the reference FASTAs and GFF annotations, and produces annotated VCFs in `snpeff/`.

### 3. Complete the analysis notebook (Steps 5–6)
Return to `Multi_species_analysis_14sp.ipynb` and run the remaining cells. This will:
- Classify mutations as synonymous, non-synonymous, or intergenic using snpEff annotations
- Perform trinucleotide context analysis
- Export results to CSV files for statistical modelling in R and figure generation

### 4. Generate figures
Open `Figures.ipynb` and run all cells. Figures are saved to `figures/`.

### 5. Run the statistical model
Open `poisson_regression_model.R` in RStudio and run all lines. Results are exported to `model_results_14sp.csv`.

---

## Input data

Raw Parsnp VCF files and reference FASTA/GFF files are stored in the shared CLIMB team directory (`/home/jovyan/shared-team/2025-masters-project/people/johanna_genomes/`) and are not included in this repository. The snpEff annotation step and the analysis notebook read directly from this shared location.

---

## Acknowledgements

Raw VCF files and reference genomes were provided by the broader mutation bias project group. snpEff database setup followed the approach used across the group. Codon enumeration logic in Step 6 was adapted from another group member's implementation.
