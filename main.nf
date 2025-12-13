#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
================================================================================
    SALIVA cfDNA ANALYSIS PIPELINE
================================================================================
    A comprehensive pipeline for saliva cell-free DNA analysis supporting:
    
    1. FASTQ MODE: Full preprocessing (merge → trim → align → dedup → blacklist)
                   followed by optional downstream analysis
    
    2. BAM MODE:   Direct downstream analysis on pre-processed BAM files
    
    Downstream Analyses Available:
      --run_nuclease      : Nuclease signature analysis (cleavage, shape, damage)
      --run_fragmentomic  : Fragmentomic ratio (short vs long fragment ratio)
      --run_endmotif      : End motif analysis (4-mer frequencies)
      --run_gquad         : G-quadruplex detection
      --run_insertsize    : Insert size distribution analysis
      --run_coverage      : Genome coverage binning (1Mb windows)
      --run_microbiome    : MetaPhlAn microbiome profiling

    Usage Examples:
      # Full pipeline from FASTQs with all downstream analyses
      nextflow run main.nf --reads 'data/*.fastq.gz' --mode saliva --run_all
      
      # Preprocessing only (no downstream)
      nextflow run main.nf --reads 'data/*.fastq.gz' --mode saliva
      
      # Downstream only from existing BAMs
      nextflow run main.nf --bam_dir /path/to/bams --run_fragmentomic --run_endmotif
      
      # Single BAM downstream analysis
      nextflow run main.nf --bam sample.bam --sample_id SAMPLE01 --run_all

    Author: UCLA Wong Laboratory
================================================================================
*/

// ============================================================================
// PARAMETERS
// ============================================================================

// ----- Input Parameters -----
// FASTQ mode inputs
params.reads          = null                    // Glob pattern for FASTQ files
params.r1             = null                    // Explicit R1 file
params.r2             = null                    // Explicit R2 file
params.r3             = null                    // Explicit R3 file (alternative to R2)

// BAM mode inputs
params.bam_dir        = null                    // Directory containing BAM files
params.bam            = null                    // Single BAM file

// Common
params.sample_id      = null                    // Sample ID (auto-derived if not set)
params.mode           = 'saliva'                // 'plasma' or 'saliva' (preprocessing mode)
params.outdir         = 'results'               // Output directory
params.reference      = null                    // Reference genome FASTA

// ----- Reference/Index Parameters -----
params.bwa_index      = null                    // BWA index prefix (plasma mode)
params.bowtie2_index  = null                    // Bowtie2 index prefix (saliva mode)
params.blacklist_bed  = null                    // Blacklist regions BED file

// ----- Downstream Analysis Flags (SALIVA-SPECIFIC) -----
// Note: These downstream analyses are designed for SALIVA cfDNA samples
// Plasma-specific downstream analyses will be added in a future update
params.run_downstream = false                   // Enable downstream analysis after preprocessing
params.run_all        = false                   // Run all downstream analyses
params.run_nuclease   = false                   // Nuclease signature analysis
params.run_fragmentomic = false                 // Fragmentomic ratio analysis
params.run_endmotif   = false                   // End motif analysis
params.run_gquad      = false                   // G-quadruplex detection
params.run_insertsize = false                   // Insert size distribution
params.run_coverage   = false                   // Genome coverage binning
params.run_microbiome = false                   // MetaPhlAn microbiome profiling
params.run_homer      = false                   // HOMER gene enrichment analysis

// ----- Downstream Analysis Flags (PLASMA-SPECIFIC) -----
// Placeholder for future plasma downstream analyses
// params.run_plasma_downstream = false

// ----- Analysis-Specific Parameters -----
params.fragmentomic_min_short   = 25
params.fragmentomic_max_short   = 100
params.fragmentomic_max_long    = 250
params.fragmentomic_bin_size    = 1000000
params.coverage_bin_size = 1000000
params.insertsize_max    = 300
params.nuclease_chunk_size  = 50000
params.nuclease_fast_mode   = false

// MetaPhlAn microbiome parameters
params.metaphlan_db      = null                 // MetaPhlAn database directory (required for --run_microbiome)
params.metaphlan_index   = null                 // MetaPhlAn database index version (optional, auto-detected)

// HOMER gene enrichment parameters
params.homer_annotations = null                 // Directory containing gene annotation BED files (required for --run_homer)
params.homer_genome      = 'hg38'               // HOMER genome (hg38, hg19, mm10, etc.)
params.homer_motif       = false                // Run findMotifsGenome (slow, optional)

// ----- Resource Parameters -----
params.qualimap_mem   = '15G'
params.max_memory     = '64.GB'
params.max_cpus       = 8
params.max_time       = '24.h'

// ============================================================================
// INCLUDE PREPROCESSING MODULES (both plasma and saliva)
// ============================================================================

include { BBMERGE        } from './modules/merge'
include { FASTP_TRIM     } from './modules/trim'
include { ALIGN_BWA      } from './modules/bwa'
include { ALIGN_BT2      } from './modules/bt2'
include { SAM2BAM        } from './modules/sam2bam'
include { QC_INT         } from './modules/qualimap_int'
include { SRSLYUMI       } from './modules/srslyumi'
include { UMITOOLS       } from './modules/umitools'
include { PICARD         } from './modules/picard'
include { CLIPPING       } from './modules/clipping'
include { AUTOSOMAL_MITO } from './modules/autosomal_mito'
include { BLACKLIST      } from './modules/blacklist'
include { QC_FINAL       } from './modules/qualimap_final'

// ============================================================================
// INCLUDE DOWNSTREAM ANALYSIS MODULES (SALIVA-SPECIFIC)
// ============================================================================
// These modules are designed for saliva cfDNA analysis
// Plasma-specific downstream modules will be added in ./modules/downstream_plasma/

include { NUCLEASE_SIGNATURE_ANALYSIS } from './modules/downstream_saliva/nuclease_signature_analysis'
include { FRAGMENTOMIC_RATIO_ANALYSIS } from './modules/downstream_saliva/fragmentomic_ratio_analysis'
include { ENDMOTIF_ANALYSIS    } from './modules/downstream_saliva/endmotif_analysis'
include { GQUAD_ANALYSIS       } from './modules/downstream_saliva/gquad_analysis'
include { INSERTSIZE_ANALYSIS  } from './modules/downstream_saliva/insertsize_analysis'
include { COVERAGE_ANALYSIS    } from './modules/downstream_saliva/coverage_analysis'
include { MICROBIOME_ANALYSIS  } from './modules/downstream_saliva/microbiome_analysis'
include { HOMER_ANALYSIS       } from './modules/downstream_saliva/homer_analysis'
include { AGGREGATE_FRAGMENTOMIC_RATIO } from './modules/downstream_saliva/aggregate'
include { AGGREGATE_ENDMOTIF   } from './modules/downstream_saliva/aggregate'
include { AGGREGATE_COVERAGE   } from './modules/downstream_saliva/aggregate'
include { AGGREGATE_GQUAD      } from './modules/downstream_saliva/aggregate'
include { AGGREGATE_INSERTSIZE } from './modules/downstream_saliva/aggregate'
include { AGGREGATE_HOMER_ANNSTATS  } from './modules/downstream_saliva/aggregate'
include { AGGREGATE_HOMER_HISTSTATS } from './modules/downstream_saliva/aggregate'

// ============================================================================
// INCLUDE DOWNSTREAM ANALYSIS MODULES (PLASMA-SPECIFIC) - PLACEHOLDER
// ============================================================================
// Uncomment and add plasma-specific modules when available
// include { PLASMA_ANALYSIS_1 } from './modules/downstream_plasma/analysis_1'

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

def deriveSampleIdFromFastq = { p ->
    new File(p.toString()).getName()
        .replaceFirst(/(?:[._-](?:R?[123]|[123])(?:_\d{3})?)\.(?:f(?:ast)?q)\.gz$/, '')
}

def deriveSampleIdFromBam = { bam_path ->
    def name = new File(bam_path.toString()).getName()
    return name.replaceFirst(/(_sorted)?(_blacklisted)?\.bam$/, '')
}

def printHeader = {
    // Determine input mode
    def inputMode = (params.bam_dir || params.bam) ? 'BAM' : 'FASTQ'
    def downstreamEnabled = params.run_all || params.run_downstream || params.run_nuclease || 
                            params.run_fragmentomic || params.run_endmotif || params.run_gquad || 
                            params.run_insertsize || params.run_coverage || params.run_homer ||
                            params.run_microbiome
    
    log.info """
    ╔═══════════════════════════════════════════════════════════════════════════╗
    ║        SALIVA cfDNA ANALYSIS PIPELINE                                     ║
    ║        UCLA Wong Laboratory                                               ║
    ╚═══════════════════════════════════════════════════════════════════════════╝
    
    Input Mode          : ${inputMode}
    ${inputMode == 'FASTQ' ? "FASTQ Pattern       : ${params.reads ?: params.r1 ?: 'Not specified'}" : "BAM Input           : ${params.bam_dir ?: params.bam ?: 'Not specified'}"}
    Processing Mode     : ${params.mode}
    Reference Genome    : ${params.reference ?: 'Not specified'}
    Output Directory    : ${params.outdir}
    
    ───────────────────────────────────────────────────────────────────────────
    PREPROCESSING (${inputMode == 'FASTQ' ? 'ENABLED' : 'SKIPPED - BAM input'})
    ───────────────────────────────────────────────────────────────────────────
    ${inputMode == 'FASTQ' ? """  ├─ Merge reads (BBMerge)
      ├─ Trim adapters (fastp)
      ├─ Align (${params.mode == 'plasma' ? 'BWA-MEM' : 'Bowtie2'})
      ├─ UMI processing (srslyumi + umi_tools)
      ├─ Deduplication (Picard)
      ${params.mode == 'plasma' ? '├─ Remove clipped reads' : '├─ (no clipping for saliva)'}
      ├─ Autosomal/Mito split
      └─ Blacklist removal""" : '  (Preprocessing skipped - using provided BAM files)'}
    
    ───────────────────────────────────────────────────────────────────────────
    DOWNSTREAM ANALYSIS - SALIVA (${downstreamEnabled ? 'ENABLED' : 'DISABLED'})
    ───────────────────────────────────────────────────────────────────────────
      ├─ Nuclease signature  : ${params.run_all || params.run_nuclease     ? '✓' : '✗'}
      ├─ Fragmentomic ratio  : ${params.run_all || params.run_fragmentomic ? '✓' : '✗'}
      ├─ End motif analysis  : ${params.run_all || params.run_endmotif  ? '✓' : '✗'}
      ├─ G-quadruplex        : ${params.run_all || params.run_gquad     ? '✓' : '✗'}
      ├─ Insert size dist    : ${params.run_all || params.run_insertsize? '✓' : '✗'}
      ├─ Coverage binning    : ${params.run_all || params.run_coverage  ? '✓' : '✗'}
      ├─ HOMER enrichment    : ${params.run_all || params.run_homer     ? '✓' : '✗'}
      └─ Microbiome (MPA)    : ${params.run_all || params.run_microbiome? '✓' : '✗'}
    
    ───────────────────────────────────────────────────────────────────────────
    DOWNSTREAM ANALYSIS - PLASMA (NOT YET AVAILABLE)
    ───────────────────────────────────────────────────────────────────────────
      └─ Coming soon...
    
    ═══════════════════════════════════════════════════════════════════════════
    """.stripIndent()
}

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {

    // Print header
    printHeader()
    
    // Determine input mode
    def isBamMode = params.bam_dir || params.bam
    def isFastqMode = params.reads || params.r1
    
    // Validate inputs
    if (!isBamMode && !isFastqMode) {
        exit 1, """
        ERROR: No input specified. Please provide either:
          FASTQ mode: --reads 'data/*.fastq.gz' or --r1 file_R1.fastq.gz --r2 file_R2.fastq.gz
          BAM mode:   --bam_dir /path/to/bams or --bam sample.bam
        """
    }
    
    if (isBamMode && isFastqMode) {
        exit 1, """
        ERROR: Cannot specify both FASTQ and BAM inputs. Please use one mode:
          FASTQ mode: --reads or --r1/--r2
          BAM mode:   --bam_dir or --bam
        """
    }
    
    // Check if downstream is requested in BAM mode
    def downstreamRequested = params.run_all || params.run_nuclease || params.run_fragmentomic || 
                              params.run_endmotif || params.run_gquad || params.run_insertsize || 
                              params.run_coverage || params.run_homer || params.run_microbiome
    
    if (isBamMode && !downstreamRequested) {
        exit 1, """
        ERROR: BAM input provided but no downstream analysis selected.
        Please specify at least one analysis:
          --run_all, --run_nuclease, --run_fragmentomic, --run_endmotif, 
          --run_gquad, --run_insertsize, --run_coverage, --run_homer, --run_microbiome
        """
    }
    
    // Warn if plasma mode with downstream (downstream is saliva-specific)
    if (params.mode == 'plasma' && downstreamRequested && isFastqMode) {
        log.warn """
        ⚠️  WARNING: Downstream analyses are currently SALIVA-SPECIFIC.
        You are running in plasma mode with downstream analysis enabled.
        Plasma-specific downstream scripts will be added in a future update.
        
        Current options:
          1. Run preprocessing only (remove --run_* flags)
          2. Switch to --mode saliva if your samples are saliva
          3. Wait for plasma downstream modules to be added
        
        Continuing with saliva downstream analyses on plasma data...
        """
    }
    
    // Validate MetaPhlAn database when microbiome is requested
    def microbiomeRequested = params.run_all || params.run_microbiome
    if (microbiomeRequested && !params.metaphlan_db) {
        log.warn """
        ⚠️  WARNING: Microbiome analysis requested but --metaphlan_db not specified.
        MetaPhlAn will attempt to use its default database location.
        
        If this fails, please:
          1. Install the MetaPhlAn database:
             metaphlan --install --bowtie2db /path/to/metaphlan_db
          
          2. Provide the path:
             --metaphlan_db /path/to/metaphlan_db
        
        On Hoffman2, use:
             --metaphlan_db /u/home/c/choi/project-wonglab/RefGenomes/mpa/
        """
    }
    
    // Reference channel
    ch_reference = params.reference ? Channel.value(file(params.reference)) : Channel.value([])
    
    // ========================================================================
    // MODE 1: BAM INPUT - Skip preprocessing, go directly to downstream
    // ========================================================================
    
    if (isBamMode) {
        log.info "Running in BAM MODE - skipping preprocessing"
        
        // Build BAM channel
        if (params.bam) {
            def sid = params.sample_id ?: deriveSampleIdFromBam(params.bam)
            ch_bam_for_downstream = Channel.of(tuple([sample_id: sid], file(params.bam)))
        } else {
            ch_bam_for_downstream = Channel
                .fromPath("${params.bam_dir}/*{_sorted_blacklisted,_blacklisted,}.bam", type: 'file')
                .map { bam ->
                    def sid = deriveSampleIdFromBam(bam)
                    tuple([sample_id: sid], bam)
                }
        }
        
        // Run downstream analyses
        RUN_DOWNSTREAM(ch_bam_for_downstream, ch_reference)
    }
    
    // ========================================================================
    // MODE 2: FASTQ INPUT - Run preprocessing, optionally run downstream
    // ========================================================================
    
    else {
        log.info "Running in FASTQ MODE - starting preprocessing"
        
        // Build FASTQ reads channel
        if (params.r1 && (params.r2 || params.r3)) {
            def sid   = params.sample_id ?: deriveSampleIdFromFastq(params.r1)
            def mate2 = params.r3 ?: params.r2
            ch_reads = Channel.of(tuple([sample_id: sid], file(params.r1), file(mate2)))
        } else {
            def RX = ~/^(.+?)(?:[._-]?(?:R?([123])|([123]))(?:_\d{3})?)\.(?:f(?:ast)?q)\.gz$/
            ch_reads = Channel
                .fromPath(params.reads, type:'file', followLinks:true)
                .map { f ->
                    def m = (f.name =~ RX)
                    assert m.matches() : "Cannot parse read filename: ${f.name}"
                    tuple(m[0][1], (m[0][2] ?: m[0][3]) as Integer, f)
                }
                .groupTuple(by:0)
                .map { id, mates, files ->
                    def i1 = mates.indexOf(1)
                    def i2 = mates.indexOf(2)
                    def i3 = mates.indexOf(3)
                    def r1 = i1 >= 0 ? files[i1] : null
                    def r2 = i2 >= 0 ? files[i2] : null
                    def r3 = i3 >= 0 ? files[i3] : null
                    def mate2 = r2 ?: r3
                    assert r1 && mate2 : "Missing mate for sample: ${id}"
                    tuple([sample_id: id], r1, mate2)
                }
        }
        
        // Run preprocessing
        ch_preprocessed = RUN_PREPROCESSING(ch_reads)
        
        // Run downstream if requested
        if (downstreamRequested || params.run_downstream) {
            RUN_DOWNSTREAM(ch_preprocessed, ch_reference)
        }
    }
}

// ============================================================================
// PREPROCESSING SUB-WORKFLOW
// ============================================================================

workflow RUN_PREPROCESSING {
    take:
    ch_reads  // tuple(meta, r1, r2)
    
    main:
    
    // 1. Merge paired reads
    ch_merge = BBMERGE(ch_reads)
    ch_merged_fastq = ch_merge.merged
    
    // 2. Adapter/quality trimming
    ch_trim = FASTP_TRIM(ch_merged_fastq)
    
    // Extract trimmed fastq for alignment
    ch_trimmed = ch_trim.trimmed.map { meta, trimmed_fastq, html, json ->
        tuple(meta, trimmed_fastq)
    }
    
    if (params.mode == 'plasma') {
        // Plasma workflow: BWA alignment
        c1 = ALIGN_BWA(ch_trimmed)
        c1_bam = c1.bam
        
        QC_INT(c1_bam)
        
        c3 = SRSLYUMI(c1_bam)
        c4 = UMITOOLS(c3.rx_bam)
        c5 = PICARD(c4.bx_bam)
        
        // Prepare for clipping: (meta, bam, bai)
        c5_for_clip = c5.deduped.map { meta, bam, bai, metrics -> tuple(meta, bam, bai) }
        c6 = CLIPPING(c5_for_clip)
        
        // Prepare for autosomal_mito: (meta, bam, [bai])
        c6_for_split = c6.clipped.map { meta, bam, bai -> tuple(meta, bam, [bai]) }
        c7 = AUTOSOMAL_MITO(c6_for_split)
        
        c8 = BLACKLIST(c7.split_bams)
        QC_FINAL(c7.split_bams)
        
        // Output: blacklisted BAM
        ch_final_bam = c8.blacklisted.map { meta, bam, bai -> tuple(meta, bam) }
    }
    else if (params.mode == 'saliva') {
        // Saliva workflow: Bowtie2 alignment
        s1 = ALIGN_BT2(ch_trimmed)
        s2 = SAM2BAM(s1.sam)
        s2_bam = s2.bam
        
        QC_INT(s2_bam)
        
        s4 = SRSLYUMI(s2_bam)
        s5 = UMITOOLS(s4.rx_bam)
        s6 = PICARD(s5.bx_bam)
        
        // Prepare for autosomal_mito (no clipping for saliva): (meta, bam, [bai, metrics])
        s6_for_split = s6.deduped.map { meta, bam, bai, metrics -> tuple(meta, bam, [bai, metrics]) }
        s7 = AUTOSOMAL_MITO(s6_for_split)
        
        s8 = BLACKLIST(s7.split_bams)
        QC_FINAL(s7.split_bams)
        
        // Output: blacklisted BAM
        ch_final_bam = s8.blacklisted.map { meta, bam, bai -> tuple(meta, bam) }
    }
    else {
        exit 1, "ERROR: --mode must be 'plasma' or 'saliva' (got: ${params.mode})"
    }
    
    emit:
    ch_final_bam
}

// ============================================================================
// DOWNSTREAM ANALYSIS SUB-WORKFLOW (SALIVA-SPECIFIC)
// ============================================================================
// Note: These analyses are designed for saliva cfDNA samples
// Plasma-specific downstream will be added as RUN_DOWNSTREAM_PLASMA

workflow RUN_DOWNSTREAM {
    take:
    ch_bam        // tuple(meta, bam)
    ch_reference  // reference fasta (can be empty)
    
    main:
    
    log.info "Starting SALIVA downstream analyses..."
    
    // Nuclease Signature Analysis
    if (params.run_all || params.run_nuclease) {
        NUCLEASE_SIGNATURE_ANALYSIS(ch_bam, ch_reference)
    }
    
    // Fragmentomic Ratio Analysis
    if (params.run_all || params.run_fragmentomic) {
        FRAGMENTOMIC_RATIO_ANALYSIS(ch_bam)
        FRAGMENTOMIC_RATIO_ANALYSIS.out.fragmentomic_csv
            .map { meta, csv -> csv }
            .collect()
            .set { ch_fragmentomic_all }
        AGGREGATE_FRAGMENTOMIC_RATIO(ch_fragmentomic_all)
    }
    
    // End Motif Analysis
    if (params.run_all || params.run_endmotif) {
        ENDMOTIF_ANALYSIS(ch_bam)
        ENDMOTIF_ANALYSIS.out.motif_ordered
            .map { meta, csv -> csv }
            .collect()
            .set { ch_endmotif_all }
        AGGREGATE_ENDMOTIF(ch_endmotif_all)
    }
    
    // G-Quadruplex Analysis
    if (params.run_all || params.run_gquad) {
        GQUAD_ANALYSIS(ch_bam, ch_reference)
        GQUAD_ANALYSIS.out.gquad_count
            .map { meta, txt -> txt }
            .collect()
            .set { ch_gquad_counts }
        GQUAD_ANALYSIS.out.gquad_total
            .map { meta, txt -> txt }
            .collect()
            .set { ch_gquad_totals }
        AGGREGATE_GQUAD(ch_gquad_counts, ch_gquad_totals)
    }
    
    // Insert Size Distribution
    if (params.run_all || params.run_insertsize) {
        INSERTSIZE_ANALYSIS(ch_bam)
        INSERTSIZE_ANALYSIS.out.insertsize_filled
            .map { meta, txt -> txt }
            .collect()
            .set { ch_insertsize_all }
        AGGREGATE_INSERTSIZE(ch_insertsize_all)
    }
    
    // Genome Coverage Binning
    if (params.run_all || params.run_coverage) {
        COVERAGE_ANALYSIS(ch_bam)
        COVERAGE_ANALYSIS.out.coverage_csv
            .map { meta, csv -> csv }
            .collect()
            .set { ch_coverage_all }
        AGGREGATE_COVERAGE(ch_coverage_all)
    }
    
    // HOMER Gene Enrichment Analysis
    if (params.run_all || params.run_homer) {
        // Check for annotations directory
        if (params.homer_annotations) {
            ch_homer_annotations = Channel.value(file(params.homer_annotations))
            
            HOMER_ANALYSIS(ch_bam, ch_homer_annotations)
            
            // Aggregate HOMER annstats results
            HOMER_ANALYSIS.out.annstats
                .map { meta, file -> file }
                .collect()
                .set { ch_homer_annstats }
            
            HOMER_ANALYSIS.out.annstats_normalized
                .map { meta, file -> file }
                .collect()
                .set { ch_homer_normalized }
            
            HOMER_ANALYSIS.out.annstats_log2ratio
                .map { meta, file -> file }
                .collect()
                .set { ch_homer_log2ratio }
            
            AGGREGATE_HOMER_ANNSTATS(ch_homer_annstats, ch_homer_normalized, ch_homer_log2ratio)
            
            // Aggregate HOMER histstats results
            HOMER_ANALYSIS.out.histstats
                .map { meta, files -> files }
                .flatten()
                .collect()
                .set { ch_homer_histstats }
            
            AGGREGATE_HOMER_HISTSTATS(ch_homer_histstats)
        } else {
            log.warn "HOMER analysis requested but --homer_annotations not provided. Skipping..."
        }
    }
    
    // Microbiome Profiling
    if (params.run_all || params.run_microbiome) {
        MICROBIOME_ANALYSIS(ch_bam)
    }
}

// ============================================================================
// WORKFLOW COMPLETION
// ============================================================================

workflow.onComplete {
    def isBamMode = params.bam_dir || params.bam
    
    log.info """
    ═══════════════════════════════════════════════════════════════════════════
    Pipeline completed ${workflow.success ? 'successfully' : 'with errors'}!
    
    Mode        : ${isBamMode ? 'BAM (downstream only)' : 'FASTQ (preprocessing + downstream)'}
    Duration    : ${workflow.duration}
    Output dir  : ${params.outdir}
    Work dir    : ${workflow.workDir}
    
    ${workflow.success ? 'All analyses completed. Check output directory for results.' : 'Check error logs for details.'}
    ═══════════════════════════════════════════════════════════════════════════
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline execution stopped with error: ${workflow.errorMessage}"
}
