# Preprocessing Modules

These modules handle the FASTQ → BAM preprocessing pipeline.

## Required Modules

Copy these from your existing pipeline or create them according to your preprocessing workflow:

| Module File | Process | Description |
|-------------|---------|-------------|
| `merge.nf` | BBMERGE | Merge paired-end reads |
| `trim.nf` | FASTP_TRIM | Adapter trimming and quality filtering |
| `bwa.nf` | ALIGN_BWA | BWA-MEM alignment (plasma mode) |
| `bt2.nf` | ALIGN_BT2 | Bowtie2 alignment (saliva mode) |
| `sam2bam.nf` | SAM2BAM | SAM to sorted BAM conversion |
| `qualimap_int.nf` | QC_INT | Intermediate QC |
| `srslyumi.nf` | SRSLYUMI | UMI tagging |
| `umitools.nf` | UMITOOLS | UMI deduplication |
| `picard.nf` | PICARD | Mark duplicates |
| `clipping.nf` | CLIPPING | Remove clipped reads (plasma only) |
| `autosomal_mito.nf` | AUTOSOMAL_MITO | Split nuclear/mitochondrial |
| `blacklist.nf` | BLACKLIST | Remove blacklist regions |
| `qualimap_final.nf` | QC_FINAL | Final QC |

## Expected Module Interface

Each module should follow this pattern:

```nextflow
process PROCESS_NAME {
    tag "${meta.sample_id}"
    
    input:
    tuple val(meta), path(input_file)
    
    output:
    tuple val(meta), path("output_file"), emit: output_name
    
    script:
    """
    # processing commands
    """
}
```

## Note

If you only need downstream analysis (BAM input mode), these preprocessing 
modules are not required. Simply use:

```bash
nextflow run main.nf --bam_dir /path/to/bams --run_all
```
