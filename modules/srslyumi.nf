/*
================================================================================
    UMI TAGGING MODULE (Flexible, Backward Compatible)
================================================================================
    Adds UMI sequences to BAM RX tag for duplicate marking.
    
    Handles multiple UMI formats automatically:
    1. UMI already in RX tag → passthrough
    2. UMI in read name (e.g., @READ:ACGTACGT) → extract to RX tag
    3. No UMI → passthrough with warning
    
    Tools used: fgbio (preferred), samtools + awk (fallback)
    Replaces non-functional 'srslyumi tag' command
    
    INPUT:  tuple(meta, bam, bai) - same as original module
    OUTPUT: tuple(meta, bam_with_rx) - BAM with RX tags added
*/

process SRSLYUMI {
    tag "${meta.sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/umi_tagged/${meta.sample_id}", mode: 'copy', enabled: params.save_intermediates ?: false
    
    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("${meta.sample_id}_RX.bam"), emit: rx_bam
    tuple val(meta), path("umi_tagging_report.txt"), emit: report, optional: true
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== UMI Tagging: ${meta.sample_id} ===" | tee umi_tagging_report.txt
    echo "Input BAM: ${bam}" | tee -a umi_tagging_report.txt
    echo "Date: \$(date)" | tee -a umi_tagging_report.txt
    echo "" | tee -a umi_tagging_report.txt
    
    # =========================================================================
    # DETECTION FUNCTIONS
    # =========================================================================
    
    # Check if BAM already has RX tags
    has_rx_tag() {
        samtools view "\$1" 2>/dev/null | head -1000 | grep -q "RX:Z:" && return 0 || return 1
    }
    
    # Check if read names contain UMI pattern
    # Common patterns: @READ:ACGTACGT, @READ_ACGTACGT, @READ+ACGTACGT
    # Also handles Illumina format: @INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y:UMI
    has_umi_in_readname() {
        local patterns_found=\$(samtools view "\$1" 2>/dev/null | head -200 | awk '{print \$1}' | grep -cE '[:_+][ACGTN]{6,12}(\$|[^ACGTN])' || echo "0")
        [[ "\$patterns_found" -gt 10 ]] && return 0 || return 1
    }
    
    # Detect UMI delimiter in read names
    detect_umi_delimiter() {
        local readname=\$(samtools view "\$1" 2>/dev/null | head -1 | awk '{print \$1}')
        if [[ "\$readname" =~ :[ACGTN]{6,12}\$ ]]; then
            echo ":"
        elif [[ "\$readname" =~ _[ACGTN]{6,12}\$ ]]; then
            echo "_"
        elif [[ "\$readname" =~ \\+[ACGTN]{6,12}\$ ]]; then
            echo "+"
        else
            echo "unknown"
        fi
    }
    
    # =========================================================================
    # SCENARIO 1: BAM already has RX tags → passthrough
    # =========================================================================
    if has_rx_tag ${bam}; then
        echo "SCENARIO: RX tags already present" | tee -a umi_tagging_report.txt
        echo "ACTION: Passthrough (no modification needed)" | tee -a umi_tagging_report.txt
        cp ${bam} ${meta.sample_id}_RX.bam
        
        # Count RX tags for report
        RX_COUNT=\$(samtools view ${bam} | head -10000 | grep -c "RX:Z:" || echo "0")
        echo "RX tags found in first 10k reads: \${RX_COUNT}" | tee -a umi_tagging_report.txt
        echo "STATUS: SUCCESS" | tee -a umi_tagging_report.txt
        exit 0
    fi
    
    # =========================================================================
    # SCENARIO 2: UMI in read name → extract to RX tag
    # =========================================================================
    if has_umi_in_readname ${bam}; then
        DELIM=\$(detect_umi_delimiter ${bam})
        echo "SCENARIO: UMI detected in read names" | tee -a umi_tagging_report.txt
        echo "DELIMITER: \${DELIM}" | tee -a umi_tagging_report.txt
        
        # Try fgbio first (most reliable)
        if command -v fgbio &> /dev/null; then
            echo "ACTION: Using fgbio CopyUmiFromReadName" | tee -a umi_tagging_report.txt
            
            fgbio CopyUmiFromReadName \\
                --input=${bam} \\
                --output=${meta.sample_id}_RX.bam \\
                2>&1 | tee -a umi_tagging_report.txt
            
            if [[ \${PIPESTATUS[0]} -eq 0 ]]; then
                echo "STATUS: SUCCESS (fgbio)" | tee -a umi_tagging_report.txt
                exit 0
            else
                echo "WARNING: fgbio failed, trying fallback" | tee -a umi_tagging_report.txt
            fi
        fi
        
        # Fallback: use samtools + awk
        echo "ACTION: Using samtools + awk fallback" | tee -a umi_tagging_report.txt
        
        samtools view -h ${bam} | awk -v delim="\${DELIM}" '
        BEGIN {OFS="\\t"; tagged=0; total=0}
        /^@/ {print; next}
        {
            total++
            # Split read name by delimiter and get last element
            n = split(\$1, parts, delim)
            umi = parts[n]
            
            # Validate UMI (6-12 bases of ACGTN)
            if (umi ~ /^[ACGTN]{6,12}\$/) {
                # Check if RX tag already exists
                has_rx = 0
                for (i=12; i<=NF; i++) {
                    if (\$i ~ /^RX:Z:/) has_rx = 1
                }
                if (!has_rx) {
                    print \$0, "RX:Z:" umi
                    tagged++
                } else {
                    print \$0
                }
            } else {
                print \$0
            }
        }
        END {
            print "Reads processed: " total > "/dev/stderr"
            print "Reads tagged: " tagged > "/dev/stderr"
        }' 2>> umi_tagging_report.txt | samtools view -bS - > ${meta.sample_id}_RX.bam
        
        echo "STATUS: SUCCESS (awk fallback)" | tee -a umi_tagging_report.txt
        exit 0
    fi
    
    # =========================================================================
    # SCENARIO 3: No UMI detected → passthrough with warning
    # =========================================================================
    echo "SCENARIO: No UMI information detected" | tee -a umi_tagging_report.txt
    echo "  - No RX tags in BAM" | tee -a umi_tagging_report.txt
    echo "  - No UMI pattern in read names" | tee -a umi_tagging_report.txt
    echo "" | tee -a umi_tagging_report.txt
    echo "ACTION: Passthrough (UMI-based deduplication will not be possible)" | tee -a umi_tagging_report.txt
    echo "WARNING: Duplicates will be marked by position only!" | tee -a umi_tagging_report.txt
    
    # Show sample of read names for debugging
    echo "" | tee -a umi_tagging_report.txt
    echo "Sample read names (first 5):" | tee -a umi_tagging_report.txt
    samtools view ${bam} | head -5 | awk '{print "  " \$1}' | tee -a umi_tagging_report.txt
    
    cp ${bam} ${meta.sample_id}_RX.bam
    echo "" | tee -a umi_tagging_report.txt
    echo "STATUS: SUCCESS (passthrough, no UMI)" | tee -a umi_tagging_report.txt
    """
}
