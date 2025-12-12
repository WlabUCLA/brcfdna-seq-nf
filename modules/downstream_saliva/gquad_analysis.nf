/*
================================================================================
    G-QUADRUPLEX ANALYSIS MODULE
================================================================================
    Detects G-quadruplex forming sequences in cfDNA fragments.
    G-quadruplexes are secondary DNA structures associated with gene regulation
    and cancer-related genomic instability.
*/

process GQUAD_ANALYSIS {
    tag "${meta.sample_id}"
    label 'process_high'
    
    publishDir "${params.outdir}/downstream/gquad/${meta.sample_id}", mode: 'copy'
    
    input:
    tuple val(meta), path(bam)
    path reference
    
    output:
    tuple val(meta), path("${meta.sample_id}_gquad.csv"), emit: gquad_results
    tuple val(meta), path("${meta.sample_id}_gquad_count.csv"), emit: gquad_count
    tuple val(meta), path("${meta.sample_id}_gquad_total.csv"), emit: gquad_total
    tuple val(meta), path("${meta.sample_id}.bed"), emit: bed_file, optional: true
    tuple val(meta), path("${meta.sample_id}.fasta"), emit: fasta_file, optional: true
    
    script:
    def has_ref = reference.size() > 0
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Convert BAM to BED
    bedtools bamtobed -i ${bam} > ${meta.sample_id}.bed
    
    # Convert BED to FASTA (requires reference)
    if [[ "${has_ref}" == "true" && -f "${reference}" ]]; then
        bedtools getfasta \\
            -fi ${reference} \\
            -bed ${meta.sample_id}.bed \\
            -fo ${meta.sample_id}.fasta
    else
        # If no reference, extract sequences from BAM directly
        samtools fasta ${bam} > ${meta.sample_id}.fasta
    fi
    
    # Run G-quadruplex regex finder
    python3 ${projectDir}/bin/fastaRegexFinder.py \\
        --fasta ${meta.sample_id}.fasta \\
        > ${meta.sample_id}_gquad.csv || true
    
    # Count G-quadruplexes found (as CSV)
    GQUAD_COUNT=\$(grep -c "chr" ${meta.sample_id}_gquad.csv || echo "0")
    echo "sample_id,gquad_count" > ${meta.sample_id}_gquad_count.csv
    echo "${meta.sample_id},\${GQUAD_COUNT}" >> ${meta.sample_id}_gquad_count.csv
    
    # Count total reads analyzed (as CSV)
    TOTAL_COUNT=\$(grep -c ">" ${meta.sample_id}.fasta || echo "0")
    echo "sample_id,total_reads" > ${meta.sample_id}_gquad_total.csv
    echo "${meta.sample_id},\${TOTAL_COUNT}" >> ${meta.sample_id}_gquad_total.csv
    
    echo "G-quadruplex analysis complete for ${meta.sample_id}"
    """
}
