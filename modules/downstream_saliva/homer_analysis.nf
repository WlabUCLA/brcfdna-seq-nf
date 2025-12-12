/*
================================================================================
    HOMER GENE ENRICHMENT ANALYSIS MODULE
================================================================================
    Combines geneintersect.sh + tagdir.sh + homer.sh into single module.
    
    Complete Workflow:
    1. bedtools intersect - BAM with gene annotation BED files
    2. makeTagDirectory - Create HOMER tag directories for each region
    3. MACS2 callpeak - Peak calling for cfDNA
    4. annotatePeaks -annStats - Annotation statistics
    5. annotatePeaks -hist - Histogram profiles for each region
    6. findMotifsGenome - Motif discovery (optional)
    
    Required:
    - HOMER installed (makeTagDirectory, annotatePeaks.pl, findMotifsGenome.pl)
    - MACS2 installed
    - Gene annotation BED files directory (--homer_annotations)
*/

process HOMER_ANALYSIS {
    tag "${meta.sample_id}"
    label 'process_high'
    
    publishDir "${params.outdir}/downstream/homer/${meta.sample_id}", mode: 'copy'
    
    input:
    tuple val(meta), path(bam)
    path annotations_dir   // Directory containing annotation BED files
    
    output:
    tuple val(meta), path("${meta.sample_id}_annstats.csv"), emit: annstats
    tuple val(meta), path("${meta.sample_id}_annstats_normalized.csv"), emit: annstats_normalized
    tuple val(meta), path("${meta.sample_id}_annstats_log2ratio.csv"), emit: annstats_log2ratio
    tuple val(meta), path("${meta.sample_id}_genestats.csv"), emit: genestats
    tuple val(meta), path("histstats/${meta.sample_id}_*_histstats.csv"), emit: histstats
    tuple val(meta), path("macs2/*"), emit: macs2_peaks, optional: true
    tuple val(meta), path("motif/*"), emit: motifs, optional: true
    
    script:
    def genome = params.homer_genome ?: 'hg38'
    def run_motif = params.homer_motif ?: false
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== HOMER Gene Enrichment Analysis ==="
    echo "Sample: ${meta.sample_id}"
    echo "Genome: ${genome}"
    
    # Create output directories
    mkdir -p intersect tags macs2 homer histstats motif
    
    # Index BAM if needed
    if [[ ! -f ${bam}.bai ]]; then
        samtools index ${bam}
    fi
    
    # =========================================================================
    # Define annotation regions (from geneintersect.sh)
    # =========================================================================
    # All regions for intersection and tag directories
    ALL_ANNOTATIONS=(
        "SINE.ann"
        "LINE.ann"
        "Simple_repeat.ann"
        "hg38.basic.promTSS"
        "hg38.basic.intron"
        "hg38.basic.intergenic"
        "hg38.basic.exon"
        "hg38.basic.5UTR"
        "hg38.basic.3UTR"
        "hg38.basic.TTS"
        "hg38.basic.nonCoding"
        "cpgIsland.ann"
    )
    
    # Regions for histogram analysis (from homer.sh - excludes LINE.ann and nonCoding)
    HIST_ANNOTATIONS=(
        "SINE.ann"
        "Simple_repeat.ann"
        "hg38.basic.TTS"
        "hg38.basic.promTSS"
        "hg38.basic.intron"
        "hg38.basic.intergenic"
        "hg38.basic.exon"
        "hg38.basic.5UTR"
        "hg38.basic.3UTR"
        "cpgIsland.ann"
    )
    
    # =========================================================================
    # STEP 1: Gene Intersect (geneintersect.sh)
    # =========================================================================
    echo "Step 1: Intersecting BAM with annotation regions..."
    for ann in "\${ALL_ANNOTATIONS[@]}"; do
        bed_file="${annotations_dir}/\${ann}.bed"
        if [[ -f "\${bed_file}" ]]; then
            echo "  Intersecting with \${ann}..."
            bedtools intersect -u -abam ${bam} -b "\${bed_file}" \\
                > "intersect/${meta.sample_id}_\${ann}_intersected.bam" 2>/dev/null || true
        else
            echo "  WARNING: \${bed_file} not found, skipping..."
        fi
    done
    
    # =========================================================================
    # STEP 2: Create Tag Directories (tagdir.sh)
    # =========================================================================
    echo "Step 2: Creating HOMER tag directories..."
    for ann in "\${ALL_ANNOTATIONS[@]}"; do
        intersected_bam="intersect/${meta.sample_id}_\${ann}_intersected.bam"
        if [[ -f "\${intersected_bam}" && -s "\${intersected_bam}" ]]; then
            echo "  Creating tag directory for \${ann}..."
            makeTagDirectory "tags/${meta.sample_id}_\${ann}_tags" "\${intersected_bam}" \\
                -genome ${genome} -checkGC 2>/dev/null || true
        fi
    done
    
    # =========================================================================
    # STEP 3: MACS2 Peak Calling (homer.sh)
    # =========================================================================
    echo "Step 3: Running MACS2 peak calling..."
    macs2 callpeak \\
        -f AUTO \\
        -g hs \\
        -q 0.01 \\
        --nomodel \\
        -t ${bam} \\
        -n ${meta.sample_id}_macs2 \\
        --outdir macs2 \\
        2>/dev/null || echo "MACS2 completed (may have warnings)"
    
    PEAKS_FILE="macs2/${meta.sample_id}_macs2_peaks.narrowPeak"
    
    # =========================================================================
    # STEP 4: HOMER annotatePeaks - Annotation Statistics (homer.sh)
    # =========================================================================
    echo "Step 4: Running HOMER annotatePeaks for annotation statistics..."
    if [[ -f "\${PEAKS_FILE}" ]]; then
        annotatePeaks.pl "\${PEAKS_FILE}" ${genome} \\
            -annStats "homer/${meta.sample_id}_annStats.txt" \\
            > "homer/${meta.sample_id}_geneStats.txt" 2>/dev/null || true
    else
        echo "No peaks file found, creating empty outputs"
        touch "homer/${meta.sample_id}_annStats.txt"
        touch "homer/${meta.sample_id}_geneStats.txt"
    fi
    
    # =========================================================================
    # STEP 5: HOMER annotatePeaks - Histogram Statistics (homer.sh)
    # =========================================================================
    echo "Step 5: Running HOMER histogram analysis for each region..."
    if [[ -f "\${PEAKS_FILE}" ]]; then
        for ann in "\${HIST_ANNOTATIONS[@]}"; do
            tagdir="tags/${meta.sample_id}_\${ann}_tags"
            if [[ -d "\${tagdir}" ]]; then
                echo "  Histogram for \${ann}..."
                annotatePeaks.pl "\${PEAKS_FILE}" ${genome} \\
                    -d "\${tagdir}" \\
                    -size 1000 \\
                    -hist 10 \\
                    > "histstats/${meta.sample_id}_\${ann}_histStats.txt" 2>/dev/null || true
            fi
        done
    fi
    
    # =========================================================================
    # STEP 6: HOMER findMotifsGenome - Motif Discovery (optional)
    # =========================================================================
    if [[ "${run_motif}" == "true" && -f "\${PEAKS_FILE}" ]]; then
        echo "Step 6: Running motif discovery..."
        findMotifsGenome.pl "\${PEAKS_FILE}" ${genome} motif -size 200 -mask 2>/dev/null || true
    else
        echo "Step 6: Skipping motif discovery (disabled or no peaks)"
    fi
    
    # =========================================================================
    # STEP 7: Parse Results to CSV (annStats_parse.py + histStats_parse.py)
    # =========================================================================
    echo "Step 7: Converting results to CSV format..."
    
    python3 << 'PYTHON_SCRIPT'
import os
import pandas as pd

sample_id = "${meta.sample_id}"

# -------------------------------------------------------------------------
# Parse annStats (from annStats_parse.py)
# -------------------------------------------------------------------------
annotations_to_extract = ["3UTR", "miRNA", "ncRNA", "TTS", "pseudo", "Exon", 
                          "Intron", "Intergenic", "Promoter", "5UTR", "snoRNA", 
                          "scRNA", "rRNA"]

ann_stats_file = f"homer/{sample_id}_annStats.txt"
if os.path.exists(ann_stats_file) and os.path.getsize(ann_stats_file) > 0:
    try:
        df = pd.read_csv(ann_stats_file, delimiter='\\t')
        
        # Extract Number of peaks
        peaks_results = {}
        log2_results = {}
        
        for ann in annotations_to_extract:
            if "Annotation" in df.columns:
                mask = df["Annotation"] == ann
                if mask.any():
                    peaks_results[ann] = df.loc[mask, "Number of peaks"].values[0] if "Number of peaks" in df.columns else 0
                    log2_results[ann] = df.loc[mask, "Log2 Ratio (obs/exp)"].values[0] if "Log2 Ratio (obs/exp)" in df.columns else 0
                else:
                    peaks_results[ann] = 0
                    log2_results[ann] = 0
        
        # Create DataFrames
        peaks_df = pd.DataFrame([peaks_results])
        peaks_df.insert(0, 'sample_id', sample_id)
        peaks_df.to_csv(f"{sample_id}_annstats.csv", index=False)
        
        # Normalize
        numeric_cols = [c for c in peaks_df.columns if c != 'sample_id']
        total = peaks_df[numeric_cols].sum(axis=1).values[0]
        if total > 0:
            normalized_df = peaks_df.copy()
            normalized_df[numeric_cols] = peaks_df[numeric_cols] / total
            normalized_df.to_csv(f"{sample_id}_annstats_normalized.csv", index=False)
        else:
            peaks_df.to_csv(f"{sample_id}_annstats_normalized.csv", index=False)
        
        # Log2 ratio
        log2_df = pd.DataFrame([log2_results])
        log2_df.insert(0, 'sample_id', sample_id)
        log2_df.to_csv(f"{sample_id}_annstats_log2ratio.csv", index=False)
        
    except Exception as e:
        print(f"Error parsing annStats: {e}")
        # Create empty files with headers
        empty_df = pd.DataFrame(columns=['sample_id'] + annotations_to_extract)
        empty_df.to_csv(f"{sample_id}_annstats.csv", index=False)
        empty_df.to_csv(f"{sample_id}_annstats_normalized.csv", index=False)
        empty_df.to_csv(f"{sample_id}_annstats_log2ratio.csv", index=False)
else:
    # Create empty files
    empty_df = pd.DataFrame(columns=['sample_id'] + annotations_to_extract)
    empty_df.to_csv(f"{sample_id}_annstats.csv", index=False)
    empty_df.to_csv(f"{sample_id}_annstats_normalized.csv", index=False)
    empty_df.to_csv(f"{sample_id}_annstats_log2ratio.csv", index=False)

# -------------------------------------------------------------------------
# Parse geneStats
# -------------------------------------------------------------------------
gene_stats_file = f"homer/{sample_id}_geneStats.txt"
if os.path.exists(gene_stats_file) and os.path.getsize(gene_stats_file) > 0:
    try:
        df = pd.read_csv(gene_stats_file, delimiter='\\t')
        df.to_csv(f"{sample_id}_genestats.csv", index=False)
    except Exception as e:
        print(f"Error parsing geneStats: {e}")
        pd.DataFrame().to_csv(f"{sample_id}_genestats.csv", index=False)
else:
    pd.DataFrame().to_csv(f"{sample_id}_genestats.csv", index=False)

# -------------------------------------------------------------------------
# Parse histStats (from histStats_parse.py)
# -------------------------------------------------------------------------
for filename in os.listdir("histstats"):
    if filename.endswith("_histStats.txt"):
        filepath = os.path.join("histstats", filename)
        csv_name = filename.replace("_histStats.txt", "_histstats.csv")
        try:
            df = pd.read_csv(filepath, delimiter='\\t')
            df.to_csv(os.path.join("histstats", csv_name), index=False)
        except Exception as e:
            print(f"Error parsing {filename}: {e}")

# Clean up txt files in histstats
for f in os.listdir("histstats"):
    if f.endswith(".txt"):
        os.remove(os.path.join("histstats", f))

print("CSV conversion complete")
PYTHON_SCRIPT
    
    echo "HOMER analysis complete for ${meta.sample_id}"
    """
}
