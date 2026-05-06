#!/bin/bash
# ============================================================================
# snpEff annotation pipeline for 14-species bacterial mutation bias project
# Run in terminal on CLIMB JupyterHub
# snpEff version: 5.4c (conda environment: ~/snpeff_env_54)
# ============================================================================

SNPEFF_ENV=~/snpeff_env_54
SNPEFF_DATA=~/snpeff_env_54/share/snpeff-5.4.0c-0/data
CONFIG=~/snpeff_env_54/share/snpeff-5.4.0c-0/snpEff.config
JOHANNA=/home/jovyan/shared-team/2025-masters-project/people/johanna_genomes
ALISON=/home/jovyan/shared-team/2025-masters-project/people/alison

# ============================================================================
# STEP 1 — Install snpEff 5.4c
# ============================================================================

conda create --prefix ~/snpeff_env_54 -y
conda install --prefix ~/snpeff_env_54 -c bioconda snpeff=5.4.0c -y

# Verify installation
conda run --prefix ~/snpeff_env_54 snpEff -version
# Expected: SnpEff 5.4c (build 2026-02-23 07:25)

# ============================================================================
# STEP 2 — Create data directories and copy reference FASTAs and GFF files
# ============================================================================

for species in \
    bordetella_pertussis \
    campylobacter_d_jejuni \
    escherichia_coli \
    haemophilus_influenzae \
    klebsiella_pneumoniae \
    listeria_monocytogenes \
    mycobacterium_tuberculosis \
    neisseria_meningitidis \
    pseudomonas_aeruginosa \
    salmonella_typhimurium \
    staphylococcus_aureus \
    staphylococcus_epidermidis \
    streptococcus_agalactiae \
    streptococcus_pneumoniae; do

    mkdir -p $SNPEFF_DATA/$species
    cp $JOHANNA/$species/reference_fasta/*.fna $SNPEFF_DATA/$species/sequences.fa
    cp $JOHANNA/$species/genomic.gff $SNPEFF_DATA/$species/genes.gff
    echo "Done: $species"
done

# ============================================================================
# STEP 3 — Add species to snpEff config file
# ============================================================================

for species in \
    bordetella_pertussis \
    campylobacter_d_jejuni \
    escherichia_coli \
    haemophilus_influenzae \
    klebsiella_pneumoniae \
    listeria_monocytogenes \
    mycobacterium_tuberculosis \
    neisseria_meningitidis \
    pseudomonas_aeruginosa \
    salmonella_typhimurium \
    staphylococcus_aureus \
    staphylococcus_epidermidis \
    streptococcus_agalactiae \
    streptococcus_pneumoniae; do

    if ! grep -q "^${species}.genome" $CONFIG; then
        echo "${species}.genome : ${species}" >> $CONFIG
        echo "Added: $species"
    else
        echo "Already present: $species"
    fi
done

# Verify config entries were added
tail -20 $CONFIG

# ============================================================================
# STEP 4 — Build snpEff databases for all 14 species
# Note: -noCheckCds and -noCheckProtein flags are required because the NCBI
# GFF files do not include cds.fa or protein.fa validation files.
# These flags affect only the final validation step, not annotation output.
# ============================================================================

for species in \
    bordetella_pertussis \
    campylobacter_d_jejuni \
    escherichia_coli \
    haemophilus_influenzae \
    klebsiella_pneumoniae \
    listeria_monocytogenes \
    mycobacterium_tuberculosis \
    neisseria_meningitidis \
    pseudomonas_aeruginosa \
    salmonella_typhimurium \
    staphylococcus_aureus \
    staphylococcus_epidermidis \
    streptococcus_agalactiae \
    streptococcus_pneumoniae; do

    echo "Building: $species"
    conda run --prefix ~/snpeff_env_54 snpEff build -gff3 -v \
        -noCheckCds -noCheckProtein \
        -c $CONFIG \
        $species 2>&1 | tail -2
    echo "---"
done

# ============================================================================
# STEP 5 — Annotate cleaned VCFs for all 14 species
# Input VCFs: vcfs/cleaned/<tag>_foranalysis.vcf
# Output:     snpeff/<tag>_annotated.vcf
# ============================================================================

mkdir -p $ALISON/snpeff

declare -A TAGS=(
    ["bordetella_pertussis"]="bpertussis"
    ["campylobacter_d_jejuni"]="cjejuni"
    ["escherichia_coli"]="ecoli"
    ["haemophilus_influenzae"]="hinfluenzae"
    ["klebsiella_pneumoniae"]="kpneumoniae"
    ["listeria_monocytogenes"]="lmonocytogenes"
    ["mycobacterium_tuberculosis"]="mtuberculosis"
    ["neisseria_meningitidis"]="nmeningitidis"
    ["pseudomonas_aeruginosa"]="paeruginosa"
    ["salmonella_typhimurium"]="styphimurium"
    ["staphylococcus_aureus"]="saureus"
    ["staphylococcus_epidermidis"]="sepidermidis"
    ["streptococcus_agalactiae"]="sagalactiae"
    ["streptococcus_pneumoniae"]="spneumoniae"
)

for folder in "${!TAGS[@]}"; do
    tag=${TAGS[$folder]}
    conda run --prefix ~/snpeff_env_54 snpEff ann -noLog \
        -c $CONFIG \
        $folder \
        $ALISON/vcfs/cleaned/${tag}_foranalysis.vcf \
        > $ALISON/snpeff/${tag}_annotated.vcf
    echo "Annotated: $folder -> ${tag}_annotated.vcf"
done

# ============================================================================
# STEP 6 — Verify output: check SNP counts per annotated VCF
# ============================================================================

for f in $ALISON/snpeff/*_annotated.vcf; do
    echo "=== $(basename $f) ==="
    grep -v "^#" $f | wc -l
done