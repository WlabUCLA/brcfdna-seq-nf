/*
================================================================================
    AGGREGATE RESULTS MODULES
================================================================================
    Combines per-sample results into multi-sample matrices for downstream
    analysis and machine learning applications.
*/

process AGGREGATE_FRAGMENTOMIC_RATIO {
    label 'process_low'
    
    publishDir "${params.outdir}/downstream/aggregated", mode: 'copy'
    
    input:
    path(fragmentomic_files)
    
    output:
    path "fragmentomic_ratio_aggregate.csv", emit: fragmentomic_matrix
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    python3 ${projectDir}/bin/fragmentomic_ratio_aggregate.py \\
        . \\
        fragmentomic_ratio_aggregate.csv
    
    echo "Fragmentomic ratio aggregation complete"
    """
}

process AGGREGATE_ENDMOTIF {
    label 'process_low'
    
    publishDir "${params.outdir}/downstream/aggregated", mode: 'copy'
    
    input:
    path(motif_files)
    
    output:
    path "endmotif_aggregate.csv", emit: endmotif_matrix
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    python3 ${projectDir}/bin/endmotif_aggregate.py \\
        . \\
        endmotif_aggregate.csv
    
    echo "End motif aggregation complete"
    """
}

process AGGREGATE_COVERAGE {
    label 'process_low'
    
    publishDir "${params.outdir}/downstream/aggregated", mode: 'copy'
    
    input:
    path(coverage_files)
    
    output:
    path "coverage_aggregate.csv", emit: coverage_original
    path "coverage_aggregate_normalized.csv", emit: coverage_normalized
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    python3 ${projectDir}/bin/coverage_aggregate.py \\
        . \\
        coverage_aggregate.csv \\
        coverage_aggregate_normalized.csv
    
    echo "Coverage aggregation complete"
    """
}

process AGGREGATE_GQUAD {
    label 'process_low'
    
    publishDir "${params.outdir}/downstream/aggregated", mode: 'copy'
    
    input:
    path(count_files)
    path(total_files)
    
    output:
    path "gquad_aggregate_counts.csv", emit: gquad_counts
    path "gquad_aggregate_totals.csv", emit: gquad_totals
    path "gquad_aggregate_ratios.csv", emit: gquad_ratios
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    python3 ${projectDir}/bin/gquad_aggregate.py \\
        . \\
        gquad_aggregate
    
    echo "G-quadruplex aggregation complete"
    """
}

process AGGREGATE_INSERTSIZE {
    label 'process_low'
    
    publishDir "${params.outdir}/downstream/aggregated", mode: 'copy'
    
    input:
    path(insertsize_files)
    
    output:
    path "insertsize_aggregate.csv", emit: insertsize_matrix
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    python3 ${projectDir}/bin/insertsize_aggregate.py \\
        . \\
        insertsize_aggregate.csv
    
    echo "Insert size aggregation complete"
    """
}

/*
================================================================================
    HOMER AGGREGATION PROCESSES
================================================================================
*/

process AGGREGATE_HOMER_ANNSTATS {
    label 'process_low'
    
    publishDir "${params.outdir}/downstream/aggregated", mode: 'copy'
    
    input:
    path(annstats_files)
    path(normalized_files)
    path(log2ratio_files)
    
    output:
    path "homer_annstats_aggregate.csv", emit: annstats_matrix
    path "homer_annstats_normalized_aggregate.csv", emit: annstats_normalized
    path "homer_annstats_log2ratio_aggregate.csv", emit: annstats_log2ratio
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    python3 << 'PYTHON_SCRIPT'
import os
import pandas as pd

# Combine annstats files
annstats_files = sorted([f for f in os.listdir('.') if f.endswith('_annstats.csv') and not 'normalized' in f and not 'log2ratio' in f])
if annstats_files:
    dfs = [pd.read_csv(f) for f in annstats_files]
    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv('homer_annstats_aggregate.csv', index=False, float_format='%.6f')
else:
    pd.DataFrame().to_csv('homer_annstats_aggregate.csv', index=False)

# Combine normalized files
norm_files = sorted([f for f in os.listdir('.') if f.endswith('_annstats_normalized.csv')])
if norm_files:
    dfs = [pd.read_csv(f) for f in norm_files]
    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv('homer_annstats_normalized_aggregate.csv', index=False, float_format='%.6f')
else:
    pd.DataFrame().to_csv('homer_annstats_normalized_aggregate.csv', index=False)

# Combine log2ratio files
log2_files = sorted([f for f in os.listdir('.') if f.endswith('_annstats_log2ratio.csv')])
if log2_files:
    dfs = [pd.read_csv(f) for f in log2_files]
    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv('homer_annstats_log2ratio_aggregate.csv', index=False, float_format='%.6f')
else:
    pd.DataFrame().to_csv('homer_annstats_log2ratio_aggregate.csv', index=False)

print("HOMER annstats aggregation complete")
PYTHON_SCRIPT
    """
}

process AGGREGATE_HOMER_HISTSTATS {
    label 'process_low'
    
    publishDir "${params.outdir}/downstream/aggregated", mode: 'copy'
    
    input:
    path(histstats_dirs)  // Directories containing *_histstats.csv files
    
    output:
    path "homer_histstats_*_aggregate.csv", emit: histstats_matrices
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    python3 << 'PYTHON_SCRIPT'
import os
import pandas as pd
from collections import defaultdict

# Collect files by region type
region_files = defaultdict(list)

# Find all histstats CSV files
for f in os.listdir('.'):
    if f.endswith('_histstats.csv'):
        # Extract region name: {sample}_{region}_histstats.csv
        for known_region in ['SINE.ann', 'LINE.ann', 'Simple_repeat.ann', 
                            'hg38.basic.promTSS', 'hg38.basic.intron', 
                            'hg38.basic.intergenic', 'hg38.basic.exon',
                            'hg38.basic.5UTR', 'hg38.basic.3UTR', 
                            'hg38.basic.TTS', 'hg38.basic.nonCoding',
                            'cpgIsland.ann']:
            if f'_{known_region}_histstats.csv' in f:
                sample = f.replace(f'_{known_region}_histstats.csv', '')
                region_files[known_region].append((sample, f))
                break

# For each region, combine second columns from all samples
for region, files in region_files.items():
    if not files:
        continue
    
    # Create range column (-500 to 500, step 10)
    range_values = list(range(-500, 501, 10))
    combined = pd.DataFrame({'Range': range_values})
    
    for sample, filepath in sorted(files):
        try:
            df = pd.read_csv(filepath)
            if df.shape[1] >= 2:
                # Extract second column (coverage/signal)
                second_col = df.iloc[:, 1]
                # Pad or truncate to match range
                second_col = second_col.reindex(range(len(range_values)), fill_value=0).reset_index(drop=True)
                combined[sample] = second_col
        except Exception as e:
            print(f"Error processing {filepath}: {e}")
    
    # Save combined file with proper CSV formatting
    safe_region = region.replace('.', '_')
    combined.to_csv(f'homer_histstats_{safe_region}_aggregate.csv', index=False, float_format='%.6f')
    print(f"Created homer_histstats_{safe_region}_aggregate.csv with {len(files)} samples")

print("HOMER histstats aggregation complete")
PYTHON_SCRIPT
    """
}
