# BRcfDNA-Seq Nextflow Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda&logoColor=white)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

> **Broad-Range cell-free DNA sequencing analysis pipeline for plasma and saliva samples**

Originally developed by Irene Choi, Neeti Swarup, and Jordan Cheng at [W Lab, UCLA](https://github.com/WlabUCLA/BRcfDNA-Seq). This Nextflow implementation provides a reproducible, scalable, and portable version of the original pipeline.

---

## 📋 Table of Contents

- [Introduction](#-introduction)
- [Computational Requirements](#-requirements)
- [Installation](#-installation)
- [Reference File Setup](#-reference-file-setup)
- [Running the Pipeline](#-running-the-pipeline)
- [Downstream Analysis Options](#-downstream-analysis-options)
- [Understanding Your Results](#-understanding-your-results)
- [Troubleshooting](#-troubleshooting)
- [Quick Reference](#-quick-reference)
- [Getting Help](#-getting-help)
- [Citation](#-citation)

---

## 🧬 Introduction

This pipeline processes **cell-free DNA (cfDNA)** from liquid biopsy samples (e.g., plasma and saliva). It will:

1. **Take in** your raw sequencing data (FASTQ files from the sequencer)
2. **Cleans and processes** the data through multiple quality-control steps to capture high quality reads
3. **Produces** analysis-ready files
4. **Analyzes** optional downstream fragmentomic and cancer detection signatures

### Inputs and Outputs

| Input | Output |
|-------|--------|
| Raw FASTQ files from sequencer (fastq.gz) | Clean BAM files ready for analysis |
| | Quality control reports |
| | Fragment size distributions |
| | Optional: Cancer signatures, motif patterns, coverage profiles |


### Pipeline Workflow Overview

```
Input: FASTQ Files (fastq.gz)
      ↓
┌─────────────────────────────────────────┐
│  STEP 1: MERGE OVERLAPPING READS        │
│  Combines paired reads that overlap     │
└─────────────────────────────────────────┘
      ↓
┌─────────────────────────────────────────┐
│  STEP 2: QUALITY TRIMMING               │
│  Removes low-quality bases & adapters   │
└─────────────────────────────────────────┘
      ↓
┌─────────────────────────────────────────┐
│  STEP 3: ALIGN TO HUMAN GENOME          │
│  Maps your reads to the reference       │
│  • Plasma mode → uses BWA-MEM           │
│  • Saliva mode → uses Bowtie2           │
└─────────────────────────────────────────┘
      ↓
┌─────────────────────────────────────────┐
│  STEP 4: REMOVE PCR DUPLICATES          │
│  Uses UMIs to identify true molecules   │
└─────────────────────────────────────────┘
      ↓
┌─────────────────────────────────────────┐
│  STEP 5: FILTER & SEPARATE              │
│  • Removes problematic regions          │
│  • Separates nuclear from mito DNA      │
└─────────────────────────────────────────┘
      ↓
┌─────────────────────────────────────────┐
│  STEP 6: DOWNSTREAM ANALYSIS (Optional) │
│  Cancer signatures, fragment patterns,  │
│  end motifs, coverage, and more         │
└─────────────────────────────────────────┘
      ↓
Output: Pre-Processed BAM files (.bam) + Analysis Results
```

---

## 💻 Computational Requirements

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| **Operating System** | Linux or macOS | Ubuntu 20.04+ or macOS 12+ |
| **RAM (Memory)** | 16 GB | 32 GB or more |
| **Free Disk Space** | 100 GB per sample | 200 GB per sample |
| **Internet Connection** | Required for setup | Required for setup |
---

## 🔧 Installation

### Option A: Using Conda (Recommended)
---

#### Step 1: Open Your Terminal

**On Mac:**
- Press `Cmd + Space`, type "Terminal", press Enter

**On Linux:**
- Press `Ctrl + Alt + T`, or find Terminal in your applications menu

---

#### Step 2: Install Miniconda (Skip if you already have conda)

Copy and paste these commands **one at a time**, pressing Enter after each:

```bash
# Download the Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer
bash Miniconda3-latest-Linux-x86_64.sh
```

**What you'll see:**
- Press **Enter** to scroll through the license agreement
- Type **yes** and press Enter to accept
- Press **Enter** to accept the default installation location
- Type **yes** when asked to initialize conda

⚠️ **Important:** After installation finishes, **close your terminal and open a new one** for the changes to take effect.

---

#### Step 3: Install Nextflow

```bash
# Download Nextflow
curl -s https://get.nextflow.io | bash

# Move it to a location where your computer can find it
sudo mv nextflow /usr/local/bin/

# Test that it works (should display a version number)
nextflow -version
```

**Expected output:** Something like `nextflow version 23.10.0`

---

#### Step 4: Download the Pipeline

```bash
# Download the pipeline code from GitHub
git clone https://github.com/WlabUCLA/BRcfDNA-Seq-nextflow.git

# Move into the pipeline folder
cd BRcfDNA-Seq-nextflow
```

---

#### Step 5: Create the Software Environment

This step installs all the bioinformatics tools the pipeline needs:

```bash
# Create the environment (this takes 10-20 minutes)
conda env create -f environment.yml

# Activate the environment
conda activate brcfdna-nf
```

**🎉 Installation complete!**

---

#### Step 6: Verify Everything Works

Run these commands to make sure everything installed correctly:

```bash
nextflow -version      # Should show: nextflow version XX.XX.X
bwa                    # Should show: Program: bwa
bowtie2 --version      # Should show: bowtie2-align-s version X.X.X
samtools --version     # Should show: samtools X.X
```

If any command shows "command not found", make sure you activated the environment with `conda activate brcfdna-nf`.

---

## 📚 Reference File Setup

The pipeline needs reference files to compare your DNA against. This is a **one-time setup**.

### Step 1: Create a Folder for References

```bash
# Create the folder structure
mkdir -p ~/references/GRCh38
cd ~/references/GRCh38
```

### Step 2: Download the Human Genome

```bash
# Download from Ensembl (this is ~800MB, takes 5-15 minutes)
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip the file
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Rename for convenience
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38.fa
```

☕ **Good time for a coffee break!**

### Step 3: Build the Alignment Index

Choose based on your sample type:

**For PLASMA samples (build BWA index):**
```bash
# This takes 1-2 hours and needs ~10GB RAM
bwa index -p GRCh38 GRCh38.fa
```

**For SALIVA samples (build Bowtie2 index):**
```bash
# This takes ~1 hour and needs ~8GB RAM
bowtie2-build GRCh38.fa GRCh38
```

💡 **Tip:** Build both indexes if you'll analyze both sample types.

### Step 4: Download the Blacklist File

The blacklist contains genomic regions that cause sequencing artifacts:

```bash
cd ~/references
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
mv hg38-blacklist.v2.bed.gz blacklist.bed
```

### Your Reference Folder Should Look Like This:

```
~/references/
├── GRCh38/
│   ├── GRCh38.fa              ← Human genome sequence
│   ├── GRCh38.amb             ← BWA index files
│   ├── GRCh38.ann
│   ├── GRCh38.bwt
│   ├── GRCh38.pac
│   ├── GRCh38.sa
│   ├── GRCh38.1.bt2           ← Bowtie2 index files
│   ├── GRCh38.2.bt2
│   ├── GRCh38.3.bt2
│   ├── GRCh38.4.bt2
│   ├── GRCh38.rev.1.bt2
│   └── GRCh38.rev.2.bt2
└── blacklist.bed      ← Blacklist regions
```

---

## 🚀 Running the Pipeline

### Before You Begin

1. Make sure you're in the pipeline folder: `cd ~/BRcfDNA-Seq-nextflow`
2. Make sure the environment is activated: `conda activate brcfdna-nf`
3. Know where your FASTQ files are located

### Understanding Your Input Files

Your sequencing facility provided FASTQ files, usually in pairs:
- `Sample_R1.fastq.gz` — Read 1
- `Sample_R2.fastq.gz` (or `Sample_R3.fastq.gz`) — Read 2

### Example 1: Analyze Plasma Samples

```bash
nextflow run main.nf \
    --reads 'data/*_R{1,3}.fastq.gz' \
    --mode plasma \
    --bwa_index ~/references/GRCh38/GRCh38 \
    --blacklist_bed ~/references/blacklist.bed \
    --outdir results/plasma_analysis \
    -profile conda
```

### Example 2: Analyze Saliva Samples

```bash
nextflow run main.nf \
    --reads 'data/*_R{1,2}.fastq.gz' \
    --mode saliva \
    --bowtie2_index ~/references/GRCh38/GRCh38 \
    --blacklist_bed ~/references/blacklist.bed \
    --outdir results/saliva_analysis \
    -profile conda
```

### Example 3: Analyze a Single Sample

```bash
nextflow run main.nf \
    --r1 data/patient001_R1.fastq.gz \
    --r2 data/patient001_R2.fastq.gz \
    --sample_id patient001 \
    --mode saliva \
    --bowtie2_index ~/references/GRCh38/GRCh38 \
    --blacklist_bed ~/references/blacklist.bed \
    --outdir results/patient001 \
    -profile conda
```

### What to Expect

- **Progress bars** will show you what's happening
- **Processing time:** 2-6 hours per sample (depending on file size)
- **Output:** A `results/` folder with all your processed files
- **Completion message** when finished

---

## 🔬 Downstream Analysis Options

After preprocessing, you can run specialized analyses on your clean BAM files.

### Available Analyses

| Analysis | Flag | What It Does |
|----------|------|--------------|
| **Nuclease Signature** | `--run_nuclease` | DNA cleavage patterns from different nucleases |
| **Fragmentomic Ratio** | `--run_fragmentomic` | Short vs long fragment ratios (cancer signature) |
| **End Motif** | `--run_endmotif` | 4-mer patterns at fragment ends |
| **G-Quadruplex** | `--run_gquad` | G4-forming sequence detection |
| **Insert Size** | `--run_insertsize` | Fragment length distribution |
| **Coverage** | `--run_coverage` | Genome-wide read coverage |
| **HOMER Enrichment** | `--run_homer` | Gene region annotation analysis |
| **Microbiome** | `--run_microbiome` | Microbial DNA profiling |
| **Run Everything** | `--run_all` | All analyses above |

### Example: Full Analysis Pipeline

```bash
nextflow run main.nf \
    --reads 'data/*_R{1,2}.fastq.gz' \
    --mode saliva \
    --bowtie2_index ~/references/GRCh38/GRCh38 \
    --blacklist_bed ~/references/blacklist.bed \
    --run_all \
    --outdir results/complete_analysis \
    -profile conda
```

### Example: Analyze Pre-existing BAM Files

If you already have processed BAM files:

```bash
nextflow run main.nf \
    --bam_dir /path/to/your/bam/files \
    --run_fragmentomic \
    --run_endmotif \
    --run_insertsize \
    -profile conda
```

### Example: HOMER Gene Enrichment

```bash
nextflow run main.nf \
    --bam_dir /path/to/bams \
    --run_homer \
    --homer_annotations /path/to/hg38_geneenrichment \
    --homer_genome hg38 \
    -profile conda
```

---

## 📊 Understanding Your Results

### Output Folder Structure

```
results/
├── merged/                          # Merged read files
├── trimmed/                         # Quality-trimmed reads
├── aligned/                         # Aligned BAM files
├── deduped/                         # After duplicate removal
├── split/                           # Separated uclear/mitochondrial reads
├── final/                           # ⭐ Final pre-processed BAM files
├── qc/                              # ⭐ Quality control reports
├── downstream/                      # Analysis results
│   ├── nuclease_signature/          # Cleavage pattern analysis
│   ├── fragmentomic_ratio/          # Fragment ratio data
│   ├── endmotif/                    # End motif frequencies
│   ├── gquad/                       # G-quadruplex results
│   ├── insertsize/                  # Fragment size data
│   ├── coverage/                    # Coverage profiles
│   ├── homer/                       # Gene enrichment results
│   ├── microbiome/                  # Microbial abundance
│   └── aggregated/                  # ⭐ Combined multi-sample data
│       ├── fragmentomic_ratio_aggregate.csv
│       ├── endmotif_aggregate.csv
│       ├── insertsize_aggregate.csv
│       ├── coverage_aggregate.csv
│       └── gquad_aggregate_ratios.csv
└── pipeline_info/                   # Execution reports & logs
```

### Important Files

| File | What It Is | Where to Find It |
|------|-----------|------------------|
| **Final BAM** | Your clean, analysis-ready data | `final/SAMPLE_sorted_blacklisted.bam` |
| **QC Report** | Interactive quality report | `qc/SAMPLE_final.QC/qualimapReport.html` |
| **Combined Results** | All samples in one table | `downstream/aggregated/*.csv` |

### How to Check Data Quality

Open the QC report in your web browser:

```bash
# On Mac
open results/qc/SAMPLE_final.QC/qualimapReport.html

# On Linux
xdg-open results/qc/SAMPLE_final.QC/qualimapReport.html
```

**What to look for:**

| Metric | ✅ Good | ⚠️ Concerning |
|--------|---------|--------------|
| **Mapping rate** | >85% (plasma), >75% (saliva) | <70% |
| **Duplication rate** | 30-60% | >70% |
| **Mean insert size** | 150-200 bp | <100 or >300 bp |

---

## 🔧 Troubleshooting

### "Command not found"

**Problem:** `nextflow: command not found`

**Solution:**
```bash
# Make sure environment is active
conda activate brcfdna-nf

# If nextflow still not found, reinstall
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### "Out of memory" Errors

**Problem:** Pipeline crashes with memory error

**Solution:**
```bash
# Run with reduced resources
nextflow run main.nf [your options] \
    --max_memory 16.GB \
    --max_cpus 4
```

### "No such file or directory"

**Problem:** Can't find reference files

**Solution:**
1. Check the path is correct: `ls ~/references/GRCh38/`
2. Use full paths starting with `/home/` or `~/`
3. Make sure files are unzipped (no `.gz` extension)

### Pipeline Stopped Halfway

**Solution:**
```bash
# Resume from where it stopped (don't start over!)
nextflow run main.nf [same options] -resume
```

### "Permission denied"

**Solution:**
```bash
chmod +x nextflow
```

---

## 📋 Quick Reference

### Essential Commands

```bash
# Activate environment (do this first!)
conda activate brcfdna-nf

# Basic plasma analysis
nextflow run main.nf \
    --reads 'data/*_R{1,3}.fastq.gz' \
    --mode plasma \
    --bwa_index ~/references/GRCh38/GRCh38 \
    --blacklist_bed ~/references/blacklist.bed \
    -profile conda

# Basic saliva analysis
nextflow run main.nf \
    --reads 'data/*_R{1,2}.fastq.gz' \
    --mode saliva \
    --bowtie2_index ~/references/GRCh38/GRCh38 \
    --blacklist_bed ~/references/blacklist.bed \
    -profile conda

# Resume a failed run
nextflow run main.nf [same options] -resume

# Clean up work files after successful run
rm -rf work/
```

### Common Parameters

| Parameter | Description |
|-----------|-------------|
| `--reads` | Path to your FASTQ files |
| `--mode` | `plasma` or `saliva` |
| `--bwa_index` | BWA index path (plasma mode) |
| `--bowtie2_index` | Bowtie2 index path (saliva mode) |
| `--blacklist_bed` | Blacklist BED file path |
| `--outdir` | Output folder name |
| `--run_all` | Run all downstream analyses |
| `-profile conda` | Use conda environment |
| `-resume` | Continue from last checkpoint |

### Downstream Analysis Flags

| Flag | Analysis |
|------|----------|
| `--run_nuclease` | Nuclease signature |
| `--run_fragmentomic` | Fragmentomic ratio |
| `--run_endmotif` | End motif analysis |
| `--run_gquad` | G-quadruplex detection |
| `--run_insertsize` | Insert size distribution |
| `--run_coverage` | Coverage analysis |
| `--run_homer` | HOMER gene enrichment |
| `--run_microbiome` | Microbiome profiling |

---

## 🆘 Getting Help

### Resources

| Resource | Best For |
|----------|----------|
| [GitHub Issues](https://github.com/WlabUCLA/BRcfDNA-Seq-nextflow/issues) | Bug reports |
| [GitHub Discussions](https://github.com/WlabUCLA/BRcfDNA-Seq-nextflow/discussions) | Questions & tips |
| Email: choi00@ucla.edu | Direct contact |

### When Reporting Issues, Include:

1. The exact command you ran
2. The complete error message
3. Your system info: `nextflow -version` and `conda --version`
4. Contents of `.nextflow.log`

---

## 📖 Citation

If you use this pipeline, please cite:

```
BRcfDNA-Seq Nextflow Pipeline
W Lab, UCLA
https://github.com/WlabUCLA/BRcfDNA-Seq-nextflow
```

### Key Tools Used

- **Nextflow:** Di Tommaso et al. (2017) Nature Biotechnology
- **BWA:** Li & Durbin (2009) Bioinformatics
- **Bowtie2:** Langmead & Salzberg (2012) Nature Methods
- **UMI-tools:** Smith et al. (2017) Genome Research
- **HOMER:** Heinz et al. (2010) Molecular Cell

---

## 🙏 Acknowledgments

- **Irene Choi, Dr. Neeti Swarup, Dr. Jordan Cheng** — Original pipeline development
- **W Lab at UCLA** — Research and protocol development
- All developers of the bioinformatics tools used

---

*Last updated: December 2025*
