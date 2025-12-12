#!/usr/bin/env python3

"""
Usage:
    python nuclease_signature.py -i <input_bam> [-r <reference_fasta>] [-o <output_dir>] [options]

Arguments:
    -i, --input       : Input BAM file (required)
    -r, --reference   : Reference genome FASTA (recommended)
    -o, --output      : Output directory (default: comprehensive_saliva_analysis)

Options:
    --max-fragments   : Maximum fragments to analyze (for testing)
    --chunk-size      : Processing chunk size (default: 50000)
    --chunked-processing : Use memory-efficient chunked processing for large files
    --disable-structures : Disable structure detection
    --disable-damage     : Disable damage analysis
    --disable-nuclease   : Disable nuclease prediction
    --disable-shape      : Disable DNA shape analysis (faster)
    --fast-mode          : Enable fast mode (larger chunks, no shape analysis)
    --workers            : Number of worker processes

Output:
    - nuclease_analysis_comprehensive.csv : Nuclease signatures
    - shape_analysis_summary.csv          : DNA shape statistics
    - comprehensive_analysis_summary.csv  : Overall summary
    - comprehensive_analysis_results.json : Full results in JSON

Example:
    python nuclease_signature.py -i sample.bam -r reference.fa -o results/
    python nuclease_signature.py -i large.bam --chunked-processing --fast-mode

Features:
    - 11 comprehensive nuclease signatures
    - DNA shape analysis (MGW, twist, roll, propeller, slide, shift)
    - Structure detection (G4, Z-DNA, hairpin, i-motif, triplex)
    - Damage analysis (oxidation, deamination, abasic sites)
"""

import pandas as pd
import numpy as np
import pysam
import argparse
import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field
from collections import defaultdict, Counter
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import json
import csv
from datetime import datetime
import sys
import time

# # Configure logging
# logging.basicConfig(
#     level=logging.INFO,
#     format='%(asctime)s - %(levelname)s - %(message)s',
#     handlers=[
#         logging.FileHandler('saliva_cfdna_analysis.log'),
#         logging.StreamHandler()
#     ]
# )
# logger = logging.getLogger(__name__)

# Configure placeholder logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.NullHandler())

# All other class and function definitions here (unchanged)...

# At the end of the file, after all class/function definitions:

def setup_logging(output_dir: Path):
    log_file = output_dir / "saliva_cfdna_analysis.log"
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    logging.info(f"Logging to: {log_file}")

# ============================================================================
# OPTIMIZED CONFIGURATION
# ============================================================================

CONFIG = {
    # Processing efficiency
    'chunk_size': 50000,          # Larger chunks for full file processing
    'max_workers': min(mp.cpu_count() - 1, 8),  # Limit parallelization
    'min_fragment_length': 30,
    'max_fragment_length': 600,
    'min_mapping_quality': 30,
    
    # Analysis toggles - ALL ENABLED for comprehensive analysis
    'enable_structure_detection': True,
    'enable_damage_analysis': True,
    'enable_nuclease_prediction': True,
    'enable_shape_analysis': True,  # ENABLED for comprehensive analysis
    
    # Memory management - MODIFIED FOR FULL FILE
    'max_fragments_in_memory': None,  # Remove memory limit for full file
    'progress_update_interval': 10000,  # Less frequent updates for performance
    
    # Simple fragment categorization
    'fragment_categories': {
        'saliva_cfDNA': (30, 600)
    },
    
    # Key oncogenic regions (reduced set)
    'oncogenic_regions': [
        ('chr17', 7571720, 7590863, 'TP53'),
        ('chr12', 25358180, 25403854, 'KRAS'),
        ('chr7', 55019017, 55211628, 'EGFR'),
        ('chr17', 43044295, 43170245, 'BRCA1'),
        ('chr13', 32315086, 32400266, 'BRCA2')
    ]
}

# Pre-compiled regex patterns for efficiency
COMPILED_PATTERNS = {
    'G4': re.compile(r'G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}'),
    'Z_DNA': re.compile(r'([CG]){6,}'),
    'hairpin': re.compile(r'([ATCG]{4,15})\w{3,20}([ATCG]{4,15})'),
    'deamination': re.compile(r'CG'),
    'oxidation': re.compile(r'GGG+'),
    'abasic': re.compile(r'[AG][AG]'),
}

# Comprehensive nuclease signatures
COMPREHENSIVE_NUCLEASE_SIGNATURES = {
    'DNASE1': {
        'preferred_cuts': ['CG', 'CA', 'TG', 'TA'],
        'cleavage_pattern': 'staggered_5p',
        'overhang_type': '5_prime',
        'overhang_length': (1, 2),
        'sequence_context': r'[CT][AG]',
        'saliva_enriched': True,
        'associated_pathway': 'saliva_secretion',
        'confidence_weight': 1.0
    },
    'DNASE1L3': {
        'preferred_cuts': ['CC', 'CT', 'TC', 'TT'],
        'cleavage_pattern': 'blunt',
        'overhang_type': 'blunt',
        'overhang_length': (0, 0),
        'sequence_context': r'[CT]{2,}',
        'saliva_enriched': True,
        'associated_pathway': 'saliva_secretion',
        'confidence_weight': 1.0
    },
    'ENDONUCLEASE_G': {
        'preferred_cuts': ['GG', 'GA', 'AG', 'AA'],
        'cleavage_pattern': 'staggered_3p',
        'overhang_type': '3_prime',
        'overhang_length': (1, 3),
        'sequence_context': r'[AG]{2,}',
        'saliva_enriched': False,
        'associated_pathway': 'apoptosis_mitochondrial',
        'confidence_weight': 0.8
    },
    'CAD': {
        'preferred_cuts': ['GA', 'GT', 'CA', 'CT'],
        'cleavage_pattern': 'internucleosomal',
        'overhang_type': '5_prime',
        'overhang_length': (2, 4),
        'sequence_context': r'[GC][AT]',
        'saliva_enriched': False,
        'associated_pathway': 'apoptosis',
        'confidence_weight': 0.9
    },
    'MICROCOCCAL_NUCLEASE': {
        'preferred_cuts': ['AT', 'TA'],
        'cleavage_pattern': 'blunt',
        'overhang_type': 'blunt',
        'overhang_length': (0, 0),
        'sequence_context': r'[AT]{2,}',
        'saliva_enriched': False,
        'associated_pathway': 'bacterial',
        'confidence_weight': 0.7
    },
    'BACTERIAL_DNASE_I': {
        'preferred_cuts': ['AT', 'TA', 'GT', 'CA'],
        'cleavage_pattern': 'random',
        'overhang_type': 'mixed',
        'overhang_length': (0, 2),
        'sequence_context': r'[AT][TA]',
        'saliva_enriched': True,
        'associated_pathway': 'bacterial_secretion',
        'confidence_weight': 0.8
    },

    'STREPTOCOCCAL_DNASE': {
        'preferred_cuts': ['GA', 'AG', 'TC', 'CT'],
        'cleavage_pattern': 'staggered_5p',
        'overhang_type': '5_prime',
        'overhang_length': (1, 3),
        'sequence_context': r'[GC][AT]',
        'saliva_enriched': True,
        'associated_pathway': 'streptococcal_virulence',
        'confidence_weight': 0.9
    },

    'STAPHYLOCOCCAL_NUCLEASE': {
        'preferred_cuts': ['CA', 'TG', 'AT', 'GC'],
        'cleavage_pattern': 'endonucleolytic',
        'overhang_type': '3_prime',
        'overhang_length': (1, 2),
        'sequence_context': r'[CA][TG]',
        'saliva_enriched': True,
        'associated_pathway': 'biofilm_formation',
        'confidence_weight': 0.8
    },

    'PSEUDOMONAS_NUCLEASE': {
        'preferred_cuts': ['GG', 'CC', 'AT'],
        'cleavage_pattern': 'random',
        'overhang_type': 'blunt',
        'overhang_length': (0, 1),
        'sequence_context': r'[GC]{2,}',
        'saliva_enriched': False,
        'associated_pathway': 'opportunistic_infection',
        'confidence_weight': 0.7
    },

    'ORAL_STREPTOCOCCUS_DNASE': {
        'preferred_cuts': ['TA', 'AT', 'CG'],
        'cleavage_pattern': 'staggered_3p',
        'overhang_type': '3_prime',
        'overhang_length': (2, 4),
        'sequence_context': r'[TA][AT]',
        'saliva_enriched': True,
        'associated_pathway': 'oral_biofilm',
        'confidence_weight': 0.9
    },

    'LACTOBACILLUS_NUCLEASE': {
        'preferred_cuts': ['AA', 'TT', 'AT'],
        'cleavage_pattern': 'blunt',
        'overhang_type': 'blunt',
        'overhang_length': (0, 0),
        'sequence_context': r'[AT]{3,}',
        'saliva_enriched': True,
        'associated_pathway': 'probiotic_activity',
        'confidence_weight': 0.6
    },

    'CANDIDA_NUCLEASE': {
        'preferred_cuts': ['GC', 'CG', 'GA', 'TC'],
        'cleavage_pattern': 'endonucleolytic',
        'overhang_type': '5_prime',
        'overhang_length': (1, 3),
        'sequence_context': r'[GC][GC]',
        'saliva_enriched': True,
        'associated_pathway': 'fungal_invasion',
        'confidence_weight': 0.7
    },

    'PREVOTELLA_NUCLEASE': {
        'preferred_cuts': ['AG', 'GA', 'CT', 'TC'],
        'cleavage_pattern': 'internucleosomal',
        'overhang_type': 'mixed',
        'overhang_length': (1, 2),
        'sequence_context': r'[AG][CT]',
        'saliva_enriched': True,
        'associated_pathway': 'periodontal_disease',
        'confidence_weight': 0.8
    },

    'ACTINOMYCES_NUCLEASE': {
        'preferred_cuts': ['CC', 'GG', 'AT'],
        'cleavage_pattern': 'random',
        'overhang_type': 'blunt',
        'overhang_length': (0, 1),
        'sequence_context': r'[CG]{2,}',
        'saliva_enriched': False,
        'associated_pathway': 'actinomycosis',
        'confidence_weight': 0.6
    },

    'FUSOBACTERIUM_NUCLEASE': {
        'preferred_cuts': ['TA', 'AT', 'GC'],
        'cleavage_pattern': 'staggered_5p',
        'overhang_type': '5_prime',
        'overhang_length': (2, 3),
        'sequence_context': r'[TA][GC]',
        'saliva_enriched': True,
        'associated_pathway': 'periodontal_pathogen',
        'confidence_weight': 0.7
    },
    'GRANZYME_A': {
        'preferred_cuts': ['AG', 'GA'],
        'cleavage_pattern': 'single_strand_nick',
        'overhang_type': 'nick',
        'overhang_length': (0, 0),
        'sequence_context': r'[AG][AG]',
        'saliva_enriched': False,
        'associated_pathway': 'cytotoxic_killing',
        'confidence_weight': 0.6
    },
    'GRANZYME_B': {
        'preferred_cuts': ['GC', 'AT'],
        'cleavage_pattern': 'internucleosomal',
        'overhang_type': '5_prime',
        'overhang_length': (2, 4),
        'sequence_context': r'[GC][AT]',
        'saliva_enriched': False,
        'associated_pathway': 'cytotoxic_killing',
        'confidence_weight': 0.6
    },
    'TREX1': {
        'preferred_cuts': ['exo'],
        'cleavage_pattern': 'exonucleolytic',
        'overhang_type': '3_prime_degraded',
        'overhang_length': (0, 0),
        'sequence_context': r'.',
        'saliva_enriched': False,
        'associated_pathway': 'cytoplasmic_DNA_sensing',
        'confidence_weight': 0.5
    },
    'AIF': {
        'preferred_cuts': ['GC', 'CG', 'AT', 'TA'],
        'cleavage_pattern': 'large_fragments',
        'overhang_type': 'mixed',
        'overhang_length': (0, 4),
        'sequence_context': r'.',
        'saliva_enriched': False,
        'associated_pathway': 'apoptosis_caspase_independent',
        'confidence_weight': 0.4
    },
    'NEUTROPHIL_ELASTASE': {
        'preferred_cuts': ['histone_cleavage'],
        'cleavage_pattern': 'NET_release',
        'overhang_type': 'not_applicable',
        'overhang_length': (0, 0),
        'sequence_context': r'.',
        'saliva_enriched': False,
        'associated_pathway': 'NETosis',
        'confidence_weight': 0.3
    },
    'OXIDATIVE_CLEAVAGE': {
        'preferred_cuts': ['GG', 'GT', 'GC'],
        'cleavage_pattern': 'random',
        'overhang_type': 'mixed',
        'overhang_length': (0, 2),
        'sequence_context': r'G[GTC]',
        'saliva_enriched': False,
        'associated_pathway': 'oxidative_stress',
        'confidence_weight': 0.7
    }
}

# Comprehensive DNA shape parameters
TETRANUCLEOTIDE_MGW = {
    # High-confidence values for minor groove width (Angstroms)
    'AAAA': 4.5, 'AAAT': 4.3, 'AACA': 4.6, 'AACG': 4.9, 'AAGA': 4.8,
    'AAGT': 4.7, 'AATA': 4.4, 'AATC': 4.6, 'AATT': 4.4, 'AAAC': 4.8,
    'AAAG': 5.0, 'AATG': 5.1,
    
    # AT-rich sequences (narrow groove)
    'ATAA': 4.6, 'ATAT': 5.2, 'ATAC': 5.4, 'ATAG': 5.3, 'ATTA': 4.9,
    'ATTT': 4.8, 'ATTC': 5.0, 'ATTG': 5.1, 'ATCA': 5.2, 'ATCT': 5.3,
    'ATCC': 5.4, 'ATCG': 5.5, 'ATGA': 5.3, 'ATGT': 5.4, 'ATGC': 5.6, 'ATGG': 5.5,
    
    # TA dinucleotides (flexible regions)
    'TAAA': 5.1, 'TAAT': 4.8, 'TAAC': 5.2, 'TAAG': 5.3, 'TATA': 5.0,
    'TATT': 4.9, 'TATC': 5.1, 'TATG': 5.2, 'TACA': 5.3, 'TACT': 5.4,
    'TACC': 5.5, 'TACG': 5.6, 'TAGA': 5.2, 'TAGT': 5.3, 'TAGC': 5.7, 'TAGG': 5.6,
    
    # GC-containing sequences (wider groove)
    'GCAA': 5.5, 'GCAT': 5.6, 'GCAC': 5.7, 'GCAG': 5.8, 'GCTA': 5.6,
    'GCTT': 5.5, 'GCTC': 5.7, 'GCTG': 5.8, 'GCCA': 5.8, 'GCCT': 5.7,
    'GCCC': 5.9, 'GCCG': 5.9, 'GCGA': 5.8, 'GCGT': 5.7, 'GCGC': 5.9, 'GCGG': 5.8,
    
    # Add more comprehensive coverage
    'CGAA': 5.4, 'CGAT': 5.5, 'CGAC': 5.6, 'CGAG': 5.7, 'CGTA': 5.8,
    'CGTT': 5.6, 'CGTC': 5.7, 'CGTG': 5.8, 'CGCA': 5.7, 'CGCT': 5.6,
    'CGCC': 5.8, 'CGCG': 5.7, 'CGGA': 5.6, 'CGGT': 5.7, 'CGGC': 5.8, 'CGGG': 5.9,
    
    'GAAA': 4.8, 'GAAT': 4.6, 'GAAC': 5.0, 'GAAG': 5.2, 'GATA': 5.5,
    'GATT': 5.3, 'GATC': 5.4, 'GATG': 5.6, 'GACA': 5.3, 'GACT': 5.4,
    'GACC': 5.5, 'GACG': 5.6, 'GAGA': 5.1, 'GAGT': 5.2, 'GAGC': 5.7, 'GAGG': 5.6,
    
    'AGAA': 4.9, 'AGAT': 5.0, 'AGAC': 5.2, 'AGAG': 5.1, 'AGTA': 5.4,
    'AGTT': 5.1, 'AGTC': 5.3, 'AGTG': 5.5, 'AGCA': 5.6, 'AGCT': 5.4,
    'AGCC': 5.7, 'AGCG': 5.8, 'AGGA': 5.3, 'AGGT': 5.4, 'AGGC': 5.8, 'AGGG': 5.7,
    
    'TGAA': 5.2, 'TGAT': 5.3, 'TGAC': 5.5, 'TGAG': 5.4, 'TGTA': 5.6,
    'TGTT': 5.4, 'TGTC': 5.5, 'TGTG': 5.6, 'TGCA': 5.7, 'TGCT': 5.5,
    'TGCC': 5.8, 'TGCG': 5.7, 'TGGA': 5.5, 'TGGT': 5.6, 'TGGC': 5.8, 'TGGG': 5.9,
    
    'GTAA': 5.4, 'GTAT': 5.5, 'GTAC': 5.6, 'GTAG': 5.5, 'GTTA': 5.3,
    'GTTT': 5.2, 'GTTC': 5.4, 'GTTG': 5.5, 'GTCA': 5.6, 'GTCT': 5.5,
    'GTCC': 5.7, 'GTCG': 5.8, 'GTGA': 5.5, 'GTGT': 5.6, 'GTGC': 5.8, 'GTGG': 5.7,
    
    'CAAA': 5.1, 'CAAT': 5.0, 'CAAC': 5.2, 'CAAG': 5.3, 'CATA': 5.4,
    'CATT': 5.2, 'CATC': 5.3, 'CATG': 5.5, 'CACA': 5.3, 'CACT': 5.4,
    'CACC': 5.5, 'CACG': 5.6, 'CAGA': 5.2, 'CAGT': 5.3, 'CAGC': 5.7, 'CAGG': 5.6,
    
    'ACAA': 5.0, 'ACAT': 5.1, 'ACAC': 5.2, 'ACAG': 5.3, 'ACTA': 5.2,
    'ACTT': 5.1, 'ACTC': 5.3, 'ACTG': 5.4, 'ACCA': 5.3, 'ACCT': 5.3,
    'ACCC': 5.4, 'ACCG': 5.5, 'ACGA': 5.4, 'ACGT': 5.3, 'ACGC': 5.6, 'ACGG': 5.5,
    
    'CTAA': 5.2, 'CTAT': 5.3, 'CTAC': 5.4, 'CTAG': 5.3, 'CTTA': 5.1,
    'CTTT': 5.0, 'CTTC': 5.2, 'CTTG': 5.3, 'CTCA': 5.4, 'CTCT': 5.3,
    'CTCC': 5.5, 'CTCG': 5.6, 'CTGA': 5.3, 'CTGT': 5.4, 'CTGC': 5.7, 'CTGG': 5.6,
    
    'TCAA': 5.1, 'TCAT': 5.2, 'TCAC': 5.3, 'TCAG': 5.4, 'TCTA': 5.3,
    'TCTT': 5.2, 'TCTC': 5.4, 'TCTG': 5.5, 'TCCA': 5.4, 'TCCT': 5.3,
    'TCCC': 5.5, 'TCCG': 5.6, 'TCGA': 5.7, 'TCGT': 5.6, 'TCGC': 5.8, 'TCGG': 5.7,
    
    'GGAA': 5.3, 'GGAT': 5.4, 'GGAC': 5.5, 'GGAG': 5.4, 'GGTA': 5.5,
    'GGTT': 5.4, 'GGTC': 5.5, 'GGTG': 5.6, 'GGCA': 5.6, 'GGCT': 5.5,
    'GGCC': 5.7, 'GGCG': 5.6, 'GGGA': 5.5, 'GGGT': 5.6, 'GGGC': 5.8, 'GGGG': 5.9,
    
    'CCAA': 5.2, 'CCAT': 5.3, 'CCAC': 5.4, 'CCAG': 5.5, 'CCTA': 5.4,
    'CCTT': 5.3, 'CCTC': 5.5, 'CCTG': 5.6, 'CCCA': 5.5, 'CCCT': 5.4,
    'CCCC': 5.6, 'CCCG': 5.7, 'CCGA': 5.6, 'CCGT': 5.5, 'CCGC': 5.7, 'CCGG': 5.8,
    
    'TTAA': 4.7, 'TTAT': 4.8, 'TTAC': 5.0, 'TTAG': 5.1, 'TTTA': 4.6,
    'TTTT': 4.5, 'TTTC': 4.8, 'TTTG': 4.9, 'TTCA': 5.0, 'TTCT': 4.9,
    'TTCC': 5.1, 'TTCG': 5.2, 'TTGA': 4.9, 'TTGT': 5.0, 'TTGC': 5.3, 'TTGG': 5.2
}

TETRANUCLEOTIDE_HELIX_TWIST = {
    # Helix twist values in degrees
    'AAAA': 35.6, 'AAAT': 35.2, 'AACA': 35.4, 'AACG': 36.2, 'AAGA': 35.8,
    'AAGT': 35.6, 'AATA': 34.8, 'AATC': 35.0, 'AATT': 34.5, 'AAAC': 35.5,
    'AAAG': 35.7, 'AATG': 35.3,
    
    'ATAA': 34.6, 'ATAT': 35.2, 'ATAC': 35.2, 'ATAG': 35.4, 'ATTA': 34.7,
    'ATTT': 34.5, 'ATTC': 34.9, 'ATTG': 35.1, 'ATCA': 35.0, 'ATCT': 35.2,
    'ATCC': 35.4, 'ATCG': 35.8, 'ATGA': 35.1, 'ATGT': 35.3, 'ATGC': 35.6, 'ATGG': 35.5,
    
    'TAAA': 36.2, 'TAAT': 36.0, 'TAAC': 36.1, 'TAAG': 36.3, 'TATA': 35.8,
    'TATT': 35.6, 'TATC': 35.9, 'TATG': 36.1, 'TACA': 36.0, 'TACT': 36.2,
    'TACC': 36.4, 'TACG': 36.6, 'TAGA': 36.1, 'TAGT': 36.3, 'TAGC': 36.5, 'TAGG': 36.4,
    
    'GCAA': 37.2, 'GCAT': 37.0, 'GCAC': 37.4, 'GCAG': 37.6, 'GCTA': 37.1,
    'GCTT': 36.9, 'GCTC': 37.3, 'GCTG': 37.5, 'GCCA': 37.8, 'GCCT': 37.6,
    'GCCC': 38.0, 'GCCG': 38.2, 'GCGA': 37.7, 'GCGT': 37.5, 'GCGC': 38.1, 'GCGG': 37.9,
    
    'CGAA': 36.8, 'CGAT': 36.6, 'CGAC': 37.0, 'CGAG': 37.2, 'CGTA': 36.7,
    'CGTT': 36.5, 'CGTC': 36.9, 'CGTG': 37.1, 'CGCA': 37.0, 'CGCT': 36.8,
    'CGCC': 37.4, 'CGCG': 36.9, 'CGGA': 37.3, 'CGGT': 37.1, 'CGGC': 37.5, 'CGGG': 37.7,
    
    'GAAA': 35.9, 'GAAT': 35.7, 'GAAC': 36.0, 'GAAG': 36.2, 'GATA': 36.1,
    'GATT': 35.8, 'GATC': 36.0, 'GATG': 36.1, 'GACA': 36.1, 'GACT': 36.3,
    'GACC': 36.5, 'GACG': 36.7, 'GAGA': 36.0, 'GAGT': 36.2, 'GAGC': 36.8, 'GAGG': 36.6,
    
    'AGAA': 35.8, 'AGAT': 35.9, 'AGAC': 36.1, 'AGAG': 36.0, 'AGTA': 36.2,
    'AGTT': 35.9, 'AGTC': 36.2, 'AGTG': 36.1, 'AGCA': 36.5, 'AGCT': 36.3,
    'AGCC': 36.7, 'AGCG': 36.9, 'AGGA': 36.1, 'AGGT': 36.3, 'AGGC': 36.8, 'AGGG': 36.6,
    
    'TGAA': 35.8, 'TGAT': 35.9, 'TGAC': 36.2, 'TGAG': 36.1, 'TGTA': 36.3,
    'TGTT': 36.0, 'TGTC': 36.2, 'TGTG': 36.4, 'TGCA': 36.6, 'TGCT': 36.3,
    'TGCC': 36.8, 'TGCG': 36.7, 'TGGA': 36.2, 'TGGT': 36.4, 'TGGC': 36.8, 'TGGG': 37.0,
    
    'GTAA': 36.0, 'GTAT': 36.1, 'GTAC': 36.3, 'GTAG': 36.2, 'GTTA': 35.9,
    'GTTT': 35.7, 'GTTC': 36.0, 'GTTG': 36.2, 'GTCA': 36.4, 'GTCT': 36.2,
    'GTCC': 36.6, 'GTCG': 36.8, 'GTGA': 36.1, 'GTGT': 36.3, 'GTGC': 36.7, 'GTGG': 36.5,
    
    'CAAA': 35.7, 'CAAT': 35.5, 'CAAC': 35.8, 'CAAG': 36.0, 'CATA': 36.0,
    'CATT': 35.7, 'CATC': 35.9, 'CATG': 36.2, 'CACA': 35.9, 'CACT': 36.1,
    'CACC': 36.3, 'CACG': 36.5, 'CAGA': 35.8, 'CAGT': 36.0, 'CAGC': 36.6, 'CAGG': 36.4,
    
    'ACAA': 35.6, 'ACAT': 35.7, 'ACAC': 35.9, 'ACAG': 36.0, 'ACTA': 35.8,
    'ACTT': 35.6, 'ACTC': 36.0, 'ACTG': 36.1, 'ACCA': 36.0, 'ACCT': 35.8,
    'ACCC': 36.2, 'ACCG': 36.4, 'ACGA': 36.1, 'ACGT': 36.0, 'ACGC': 36.5, 'ACGG': 36.3,
    
    'CTAA': 35.8, 'CTAT': 35.9, 'CTAC': 36.1, 'CTAG': 36.0, 'CTTA': 35.7,
    'CTTT': 35.5, 'CTTC': 35.8, 'CTTG': 36.0, 'CTCA': 36.1, 'CTCT': 35.9,
    'CTCC': 36.3, 'CTCG': 36.5, 'CTGA': 36.0, 'CTGT': 36.1, 'CTGC': 36.6, 'CTGG': 36.4,
    
    'TCAA': 35.7, 'TCAT': 35.8, 'TCAC': 36.0, 'TCAG': 36.1, 'TCTA': 35.9,
    'TCTT': 35.7, 'TCTC': 36.1, 'TCTG': 36.2, 'TCCA': 36.1, 'TCCT': 35.9,
    'TCCC': 36.3, 'TCCG': 36.5, 'TCGA': 36.6, 'TCGT': 36.4, 'TCGC': 36.8, 'TCGG': 36.6,
    
    'GGAA': 36.1, 'GGAT': 36.2, 'GGAC': 36.4, 'GGAG': 36.3, 'GGTA': 36.4,
    'GGTT': 36.1, 'GGTC': 36.3, 'GGTG': 36.5, 'GGCA': 36.7, 'GGCT': 36.4,
    'GGCC': 36.9, 'GGCG': 36.6, 'GGGA': 36.3, 'GGGT': 36.5, 'GGGC': 36.8, 'GGGG': 37.0,
    
    'CCAA': 35.9, 'CCAT': 36.0, 'CCAC': 36.2, 'CCAG': 36.3, 'CCTA': 36.1,
    'CCTT': 35.9, 'CCTC': 36.3, 'CCTG': 36.5, 'CCCA': 36.3, 'CCCT': 36.1,
    'CCCC': 36.5, 'CCCG': 36.7, 'CCGA': 36.4, 'CCGT': 36.2, 'CCGC': 36.6, 'CCGG': 36.8,
    
    'TTAA': 34.8, 'TTAT': 34.9, 'TTAC': 35.1, 'TTAG': 35.2, 'TTTA': 34.7,
    'TTTT': 34.5, 'TTTC': 34.9, 'TTTG': 35.0, 'TTCA': 35.1, 'TTCT': 34.9,
    'TTCC': 35.3, 'TTCG': 35.5, 'TTGA': 35.0, 'TTGT': 35.1, 'TTGC': 35.6, 'TTGG': 35.4
}

TETRANUCLEOTIDE_ROLL = {
    # Roll values in degrees (capped to realistic range ±5°)
    'AAAA': 0.8, 'AAAT': 0.6, 'AACA': 0.4, 'AACG': -0.2, 'AAGA': 0.2,
    'AAGT': 0.0, 'AATA': 1.2, 'AATC': 0.8, 'AATT': 1.0, 'AAAC': 0.3,
    'AAAG': 0.1, 'AATG': 0.5,
    
    'ATAA': 2.1, 'ATAT': 1.8, 'ATAC': 1.5, 'ATAG': 1.2, 'ATTA': 2.3,
    'ATTT': 2.5, 'ATTC': 2.0, 'ATTG': 1.7, 'ATCA': 1.6, 'ATCT': 1.4,
    'ATCC': 1.1, 'ATCG': 0.8, 'ATGA': 1.8, 'ATGT': 1.5, 'ATGC': 1.2, 'ATGG': 1.0,
    
    'TAAA': 2.0, 'TAAT': 1.8, 'TAAC': 1.6, 'TAAG': 1.4, 'TATA': 2.2,
    'TATT': 2.4, 'TATC': 1.9, 'TATG': 1.6, 'TACA': 1.5, 'TACT': 1.3,
    'TACC': 1.0, 'TACG': 0.7, 'TAGA': 1.7, 'TAGT': 1.4, 'TAGC': 1.1, 'TAGG': 0.9,
    
    'GCAA': -0.5, 'GCAT': -0.3, 'GCAC': -0.7, 'GCAG': -0.9, 'GCTA': -0.4,
    'GCTT': -0.2, 'GCTC': -0.6, 'GCTG': -0.8, 'GCCA': -0.8, 'GCCT': -0.6,
    'GCCC': -1.0, 'GCCG': -1.2, 'GCGA': -0.7, 'GCGT': -0.5, 'GCGC': -1.1, 'GCGG': -0.9,
    
    'CGAA': -0.4, 'CGAT': -0.2, 'CGAC': -0.6, 'CGAG': -0.8, 'CGTA': -0.3,
    'CGTT': -0.1, 'CGTC': -0.5, 'CGTG': -0.7, 'CGCA': -0.6, 'CGCT': -0.4,
    'CGCC': -1.0, 'CGCG': -0.5, 'CGGA': -0.9, 'CGGT': -0.7, 'CGGC': -1.1, 'CGGG': -1.3,
    
    'GAAA': 0.7, 'GAAT': 0.5, 'GAAC': 0.3, 'GAAG': 0.1, 'GATA': 0.9,
    'GATT': 0.6, 'GATC': 0.4, 'GATG': 0.2, 'GACA': 0.4, 'GACT': 0.2,
    'GACC': -0.1, 'GACG': -0.3, 'GAGA': 0.6, 'GAGT': 0.4, 'GAGC': -0.2, 'GAGG': -0.4,
    
    'AGAA': 0.6, 'AGAT': 0.4, 'AGAC': 0.2, 'AGAG': 0.5, 'AGTA': 0.8,
    'AGTT': 0.5, 'AGTC': 0.3, 'AGTG': 0.1, 'AGCA': 0.1, 'AGCT': 0.3,
    'AGCC': -0.2, 'AGCG': -0.4, 'AGGA': 0.4, 'AGGT': 0.2, 'AGGC': -0.3, 'AGGG': -0.5,
    
    'TGAA': 0.8, 'TGAT': 0.6, 'TGAC': 0.4, 'TGAG': 0.7, 'TGTA': 1.0,
    'TGTT': 0.7, 'TGTC': 0.5, 'TGTG': 0.3, 'TGCA': 0.3, 'TGCT': 0.5,
    'TGCC': 0.0, 'TGCG': -0.2, 'TGGA': 0.6, 'TGGT': 0.4, 'TGGC': -0.1, 'TGGG': -0.3,
    
    'GTAA': 0.5, 'GTAT': 0.7, 'GTAC': 0.4, 'GTAG': 0.6, 'GTTA': 0.4,
    'GTTT': 0.2, 'GTTC': 0.5, 'GTTG': 0.7, 'GTCA': 0.2, 'GTCT': 0.4,
    'GTCC': -0.1, 'GTCG': -0.3, 'GTGA': 0.6, 'GTGT': 0.4, 'GTGC': -0.2, 'GTGG': -0.4,
    
    'CAAA': 0.3, 'CAAT': 0.1, 'CAAC': 0.4, 'CAAG': 0.6, 'CATA': 0.6,
    'CATT': 0.3, 'CATC': 0.5, 'CATG': 0.8, 'CACA': 0.5, 'CACT': 0.7,
    'CACC': 0.2, 'CACG': 0.0, 'CAGA': 0.1, 'CAGT': 0.3, 'CAGC': -0.1, 'CAGG': -0.3,
    
    'ACAA': 0.2, 'ACAT': 0.4, 'ACAC': 0.1, 'ACAG': 0.3, 'ACTA': 0.4,
    'ACTT': 0.1, 'ACTC': 0.5, 'ACTG': 0.7, 'ACCA': 0.1, 'ACCT': 0.3,
    'ACCC': -0.2, 'ACCG': -0.4, 'ACGA': 0.3, 'ACGT': 0.1, 'ACGC': -0.3, 'ACGG': -0.5,
    
    'CTAA': 0.4, 'CTAT': 0.6, 'CTAC': 0.3, 'CTAG': 0.5, 'CTTA': 0.2,
    'CTTT': 0.0, 'CTTC': 0.3, 'CTTG': 0.5, 'CTCA': 0.3, 'CTCT': 0.5,
    'CTCC': 0.0, 'CTCG': -0.2, 'CTGA': 0.1, 'CTGT': 0.3, 'CTGC': -0.2, 'CTGG': -0.4,
    
    'TCAA': 0.1, 'TCAT': 0.3, 'TCAC': 0.0, 'TCAG': 0.2, 'TCTA': 0.3,
    'TCTT': 0.0, 'TCTC': 0.4, 'TCTG': 0.6, 'TCCA': 0.0, 'TCCT': 0.2,
    'TCCC': -0.3, 'TCCG': -0.5, 'TCGA': 0.2, 'TCGT': 0.0, 'TCGC': -0.4, 'TCGG': -0.6,
    
    'GGAA': 0.2, 'GGAT': 0.4, 'GGAC': 0.1, 'GGAG': 0.3, 'GGTA': 0.4,
    'GGTT': 0.1, 'GGTC': 0.3, 'GGTG': 0.5, 'GGCA': 0.1, 'GGCT': 0.3,
    'GGCC': -0.2, 'GGCG': -0.4, 'GGGA': 0.3, 'GGGT': 0.1, 'GGGC': -0.3, 'GGGG': -0.5,
    
    'CCAA': 0.0, 'CCAT': 0.2, 'CCAC': -0.1, 'CCAG': 0.1, 'CCTA': 0.2,
    'CCTT': -0.1, 'CCTC': 0.3, 'CCTG': 0.5, 'CCCA': -0.1, 'CCCT': 0.1,
    'CCCC': -0.4, 'CCCG': -0.6, 'CCGA': 0.1, 'CCGT': -0.1, 'CCGC': -0.5, 'CCGG': -0.7,
    
    'TTAA': 1.5, 'TTAT': 1.7, 'TTAC': 1.2, 'TTAG': 1.4, 'TTTA': 1.3,
    'TTTT': 1.1, 'TTTC': 1.4, 'TTTG': 1.6, 'TTCA': 1.2, 'TTCT': 1.4,
    'TTCC': 0.9, 'TTCG': 0.7, 'TTGA': 1.6, 'TTGT': 1.4, 'TTGC': 1.0, 'TTGG': 0.8
}

# Propeller Twist Parameters (Dinucleotide)
PROPELLER_TWIST = {
    'AA': -12.5, 'AG': -10.8, 'GA': -11.2, 'GG': -9.8,
    'AT': -14.8, 'AC': -8.9, 'GT': -10.1, 'GC': -7.2,
    'TA': -15.2, 'TG': -11.5, 'CA': -9.7, 'CG': -6.8,
    'TT': -12.5, 'TC': -10.3, 'CT': -9.4, 'CC': -8.1
}

# Slide and Shift Parameters (Dinucleotide)
SLIDE_PARAMETERS = {
    'AA': -0.15, 'AT': -0.84, 'AG': 0.02, 'AC': 0.13,
    'TA': 0.28, 'TT': -0.15, 'TG': 0.13, 'TC': 0.02,
    'GA': -0.02, 'GT': 0.10, 'GG': -0.08, 'GC': 0.05,
    'CA': 0.18, 'CT': -0.05, 'CG': 0.09, 'CC': -0.12
}

SHIFT_PARAMETERS = {
    'AA': 0.04, 'AT': -0.03, 'AG': -0.01, 'AC': 0.18,
    'TA': 0.15, 'TT': 0.02, 'TG': -0.08, 'TC': 0.11,
    'GA': 0.07, 'GT': -0.04, 'GG': 0.06, 'GC': 0.03,
    'CA': -0.06, 'CT': 0.12, 'CG': -0.02, 'CC': 0.08
}

# ============================================================================
# OPTIMIZED FRAGMENT CLASS
# ============================================================================

@dataclass
class OptimizedSalivaFragment:
    """Enhanced fragment class with comprehensive analysis features"""
    # Essential data
    read_id: str
    chromosome: str
    start: int
    end: int
    length: int
    strand: str
    mapping_quality: int
    
    # Optional data
    sequence: Optional[str] = None
    gc_content: Optional[float] = None
    
    # Nuclease analysis (enhanced)
    predicted_nuclease: Optional[str] = None
    nuclease_confidence: float = 0.0
    predicted_5p_nuclease: Optional[str] = None
    predicted_3p_nuclease: Optional[str] = None
    
    # Structure flags
    has_g4: bool = False
    has_z_dna: bool = False
    has_hairpin: bool = False
    has_i_motif: bool = False
    has_triplex: bool = False
    
    # Shape parameters (comprehensive)
    mean_mgw: float = 0.0
    mean_helix_twist: float = 0.0
    mean_roll: float = 0.0
    mean_propeller_twist: float = 0.0
    mean_slide: float = 0.0
    mean_shift: float = 0.0
    
    # Shape-derived metrics
    bendability: float = 0.0
    flexibility_score: float = 0.0
    thermal_stability: float = 0.0
    nuclease_accessibility: float = 0.0
    
    # Damage metrics
    damage_score: float = 0.0
    oxidation_count: int = 0
    deamination_count: int = 0
    abasic_count: int = 0
    
    # Regional annotation
    in_oncogenic_region: bool = False
    oncogene: Optional[str] = None
    
    # Saliva-specific
    saliva_nuclease_enriched: bool = False
    potential_bacterial_origin: bool = False

# ============================================================================
# OPTIMIZED ANALYZERS
# ============================================================================

class ComprehensiveShapeCalculator:
    """Comprehensive DNA shape analysis with full parameter coverage"""
    
    def __init__(self):
        self.mgw_table = TETRANUCLEOTIDE_MGW
        self.twist_table = TETRANUCLEOTIDE_HELIX_TWIST
        self.roll_table = TETRANUCLEOTIDE_ROLL
        self.propeller_table = PROPELLER_TWIST
        self.slide_table = SLIDE_PARAMETERS
        self.shift_table = SHIFT_PARAMETERS
        
    def calculate_shape_profile(self, sequence: str) -> Dict[str, Any]:
        """Calculate comprehensive shape parameters for a sequence"""
        if not sequence or len(sequence) < 4:
            return self._empty_shape_profile()
        
        sequence = sequence.upper()
        
        # Calculate tetranucleotide-based parameters
        mgw_values = []
        twist_values = []
        roll_values = []
        
        for i in range(len(sequence) - 3):
            tetra = sequence[i:i+4]
            if 'N' not in tetra:
                mgw_values.append(self.mgw_table.get(tetra, 5.5))
                twist_values.append(self.twist_table.get(tetra, 36.0))
                roll_values.append(self.roll_table.get(tetra, 0.0))
        
        # Calculate dinucleotide-based parameters
        propeller_values = []
        slide_values = []
        shift_values = []
        
        for i in range(len(sequence) - 1):
            dimer = sequence[i:i+2]
            if 'N' not in dimer:
                propeller_values.append(self.propeller_table.get(dimer, -10.0))
                slide_values.append(self.slide_table.get(dimer, 0.0))
                shift_values.append(self.shift_table.get(dimer, 0.0))
        
        # Calculate derived metrics
        bendability = self._calculate_bendability(sequence, roll_values, mgw_values)
        flexibility = self._calculate_flexibility(sequence, twist_values, propeller_values)
        thermal_stability = self._calculate_thermal_stability(sequence, mgw_values, twist_values)
        nuclease_accessibility = self._calculate_nuclease_accessibility(mgw_values, bendability, flexibility)
        
        return {
            'mean_mgw': np.mean(mgw_values) if mgw_values else 0.0,
            'mean_helix_twist': np.mean(twist_values) if twist_values else 0.0,
            'mean_roll': np.mean(roll_values) if roll_values else 0.0,
            'mean_propeller_twist': np.mean(propeller_values) if propeller_values else 0.0,
            'mean_slide': np.mean(slide_values) if slide_values else 0.0,
            'mean_shift': np.mean(shift_values) if shift_values else 0.0,
            'bendability': bendability,
            'flexibility_score': flexibility,
            'thermal_stability': thermal_stability,
            'nuclease_accessibility': nuclease_accessibility
        }
    
    def _calculate_bendability(self, sequence: str, roll_values: List[float], mgw_values: List[float]) -> float:
        """Calculate DNA bendability from roll and MGW"""
        if not roll_values or not mgw_values:
            return 0.0
        
        # Roll variance contribution
        roll_variance = np.var(roll_values) if len(roll_values) > 1 else 0
        roll_score = min(roll_variance / 2.0, 1.0)
        
        # MGW contribution (narrow groove = more bendable)
        mgw_mean = np.mean(mgw_values)
        mgw_score = max(0, (6.0 - mgw_mean) / 2.0)
        
        # A-tract contribution
        a_tract_score = len(re.findall(r'A{4,}', sequence)) / (len(sequence) / 10)
        a_tract_score = min(a_tract_score, 1.0)
        
        # TATA box contribution
        tata_score = 1.0 if 'TATA' in sequence else 0.0
        
        return min(roll_score * 0.4 + mgw_score * 0.3 + a_tract_score * 0.2 + tata_score * 0.1, 1.0)
    
    def _calculate_flexibility(self, sequence: str, twist_values: List[float], propeller_values: List[float]) -> float:
        """Calculate DNA flexibility"""
        if not twist_values:
            return 0.0
        
        # Twist variance
        twist_var = np.var(twist_values) if len(twist_values) > 1 else 0
        twist_flex = min(twist_var / 10.0, 1.0)
        
        # TA/TG content (flexible steps)
        flexible_steps = sequence.count('TA') + sequence.count('TG')
        step_flex = min(flexible_steps / (len(sequence) / 4), 1.0)
        
        # Propeller twist variance
        prop_var = np.var(propeller_values) if len(propeller_values) > 1 else 0
        prop_flex = min(prop_var / 20.0, 1.0)
        
        return twist_flex * 0.4 + step_flex * 0.3 + prop_flex * 0.3
    
    def _calculate_thermal_stability(self, sequence: str, mgw_values: List[float], twist_values: List[float]) -> float:
        """Calculate thermal stability"""
        # GC content contribution
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
        
        # Base stacking contribution
        strong_stacks = ['GC', 'CG', 'GG', 'CC']
        stack_score = sum(1 for i in range(len(sequence)-1) if sequence[i:i+2] in strong_stacks)
        stack_score = stack_score / (len(sequence) - 1) if len(sequence) > 1 else 0
        
        # Shape contribution
        shape_score = 0.5
        if mgw_values and twist_values:
            # More regular structure = more stable
            mgw_regularity = 1.0 - (np.std(mgw_values) / np.mean(mgw_values)) if np.mean(mgw_values) > 0 else 0
            twist_regularity = 1.0 - (np.std(twist_values) / np.mean(twist_values)) if np.mean(twist_values) > 0 else 0
            shape_score = (mgw_regularity + twist_regularity) / 2
        
        return min(gc_content * 0.4 + stack_score * 0.3 + shape_score * 0.3, 1.0)
    
    def _calculate_nuclease_accessibility(self, mgw_values: List[float], bendability: float, flexibility: float) -> float:
        """Calculate nuclease accessibility"""
        if not mgw_values:
            return 0.5
        
        # Wider groove = more accessible
        mgw_score = np.mean(mgw_values) / 6.0
        
        # Higher bendability and flexibility = more accessible
        combined_score = (mgw_score * 0.4 + bendability * 0.3 + flexibility * 0.3)
        
        return min(max(combined_score, 0.0), 1.0)
    
    def _empty_shape_profile(self) -> Dict[str, Any]:
        """Return empty shape profile"""
        return {
            'mean_mgw': 0.0,
            'mean_helix_twist': 0.0,
            'mean_roll': 0.0,
            'mean_propeller_twist': 0.0,
            'mean_slide': 0.0,
            'mean_shift': 0.0,
            'bendability': 0.0,
            'flexibility_score': 0.0,
            'thermal_stability': 0.0,
            'nuclease_accessibility': 0.0
        }

class EnhancedStructureDetector:
    """Enhanced structure detection with more comprehensive patterns"""
    
    def __init__(self):
        self.patterns = COMPILED_PATTERNS
        # Add more comprehensive patterns
        self.patterns['i_motif'] = re.compile(r'C{3,}\w{1,7}C{3,}\w{1,7}C{3,}\w{1,7}C{3,}')
        self.patterns['triplex'] = re.compile(r'[AG]{8,}|[CT]{8,}')
    
    def detect_structures(self, sequence: str) -> Dict[str, bool]:
        """Enhanced structure detection"""
        if not sequence or len(sequence) < 10:
            return {'has_g4': False, 'has_z_dna': False, 'has_hairpin': False, 
                   'has_i_motif': False, 'has_triplex': False}
        
        sequence = sequence.upper()
        
        return {
            'has_g4': bool(self.patterns['G4'].search(sequence)),
            'has_z_dna': bool(self.patterns['Z_DNA'].search(sequence)),
            'has_hairpin': bool(self.patterns['hairpin'].search(sequence)),
            'has_i_motif': bool(self.patterns['i_motif'].search(sequence)),
            'has_triplex': bool(self.patterns['triplex'].search(sequence))
        }

class EnhancedDamageAnalyzer:
    """Enhanced damage analysis with more damage types"""
    
    def __init__(self):
        self.patterns = COMPILED_PATTERNS
        # Add more damage patterns
        self.patterns['strand_breaks'] = re.compile(r'[AT]{4,}')
        self.patterns['uv_damage'] = re.compile(r'TT|TC')
    
    def analyze_damage(self, sequence: str) -> Dict[str, Any]:
        """Enhanced damage assessment"""
        if not sequence:
            return {'damage_score': 0.0, 'oxidation_count': 0, 'deamination_count': 0, 'abasic_count': 0}
        
        sequence = sequence.upper()
        
        # Count various damage types
        oxidation_count = len(self.patterns['oxidation'].findall(sequence))
        deamination_count = len(self.patterns['deamination'].findall(sequence))
        abasic_count = len(re.findall(r'[AG][AG]', sequence))  # Purine runs prone to depurination
        uv_damage = len(self.patterns['uv_damage'].findall(sequence))
        
        # Enhanced damage score
        total_damage = oxidation_count + deamination_count + abasic_count + uv_damage
        damage_score = total_damage / len(sequence)
        
        return {
            'damage_score': damage_score,
            'oxidation_count': oxidation_count,
            'deamination_count': deamination_count,
            'abasic_count': abasic_count,
            'uv_damage_count': uv_damage
        }

class ComprehensiveNucleasePredictor:
    """Comprehensive nuclease prediction with full signature database"""
    
    def __init__(self):
        self.signatures = COMPREHENSIVE_NUCLEASE_SIGNATURES
        # Pre-compile regex patterns for efficiency
        self.compiled_patterns = {}
        for nuclease, sig in self.signatures.items():
            try:
                self.compiled_patterns[nuclease] = re.compile(sig['sequence_context'])
            except:
                self.compiled_patterns[nuclease] = None
    
    def predict_nuclease(self, sequence: str, overhang_5p: int = 0, overhang_3p: int = 0) -> Tuple[Optional[str], float, Optional[str], float]:
        """Predict nucleases for both 5' and 3' ends"""
        if not sequence or len(sequence) < 4:
            return None, 0.0, None, 0.0
        
        # Get end motifs
        motif_5p = sequence[:2].upper()
        motif_3p = sequence[-2:].upper()
        
        # Predict for each end
        nuclease_5p, conf_5p = self._predict_single_end(motif_5p, overhang_5p, sequence)
        nuclease_3p, conf_3p = self._predict_single_end(motif_3p, overhang_3p, sequence)
        
        return nuclease_5p, conf_5p, nuclease_3p, conf_3p
    
    def _predict_single_end(self, motif: str, overhang: int, full_sequence: str) -> Tuple[Optional[str], float]:
        """Predict nuclease for single end"""
        scores = {}
        
        for nuclease, sig in self.signatures.items():
            score = 0.0
            
            # Motif preference
            if motif in sig['preferred_cuts']:
                score += 0.4 * sig['confidence_weight']
            
            # Overhang pattern matching
            if sig['overhang_type'] == 'blunt' and overhang == 0:
                score += 0.3
            elif sig['overhang_type'] in ['5_prime', '3_prime'] and sig['overhang_length'][0] <= overhang <= sig['overhang_length'][1]:
                score += 0.3
            
            # Sequence context
            if self.compiled_patterns[nuclease] and self.compiled_patterns[nuclease].search(full_sequence):
                score += 0.2
            
            # Saliva enrichment bonus
            if sig['saliva_enriched']:
                score += 0.1
            
            scores[nuclease] = score
        
        # Get best prediction
        if scores:
            best_nuclease, best_score = max(scores.items(), key=lambda x: x[1])
            if best_score > 0.3:  # Minimum confidence threshold
                return best_nuclease, best_score
        
        return None, 0.0

# ============================================================================
# OPTIMIZED FRAGMENT PROCESSOR
# ============================================================================

class OptimizedFragmentProcessor:
    """High-performance fragment processor"""
    
    def __init__(self, reference_fasta: Optional[str] = None):
        self.reference = None
        if reference_fasta:
            try:
                self.reference = pysam.FastaFile(reference_fasta)
                logger.info("Reference genome loaded successfully")
            except Exception as e:
                logger.warning(f"Could not load reference: {e}")
        

        # Initialize analyzers
        self.structure_detector = EnhancedStructureDetector() if CONFIG['enable_structure_detection'] else None
        self.damage_analyzer = EnhancedDamageAnalyzer() if CONFIG['enable_damage_analysis'] else None
        self.nuclease_predictor = ComprehensiveNucleasePredictor() if CONFIG['enable_nuclease_prediction'] else None



        # Pre-compute oncogenic regions for fast lookup
        self.oncogenic_lookup = self._build_oncogenic_lookup()
        
        logger.info("Optimized processor initialized")
    
    def _build_oncogenic_lookup(self) -> Dict:
        """Build fast lookup for oncogenic regions"""
        lookup = defaultdict(list)
        for chrom, start, end, gene in CONFIG['oncogenic_regions']:
            lookup[chrom].append((start, end, gene))
        return dict(lookup)
    
    def process_bam_file_chunked(self, bam_path: str, output_dir: Path, 
                                chunk_size: int = 100000) -> Dict[str, Any]:
        """Process BAM file in chunks for very large files"""
        logger.info(f"Processing large BAM file in chunks of {chunk_size}")
        
        # Initialize summary statistics
        total_processed = 0
        total_skipped = 0
        chunk_count = 0
        
        # Prepare CSV output file
        fragment_csv = output_dir / 'saliva_fragments_chunked.csv'
        csv_header_written = False
        
        try:
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                chunk_fragments = []
                
                progress_bar = tqdm(
                    bam.fetch(), 
                    desc="Processing reads (chunked)",
                    unit="reads"
                )
                
                for read in progress_bar:
                    # Quick filters
                    if (read.is_unmapped or 
                        read.mapping_quality < CONFIG['min_mapping_quality'] or
                        read.is_duplicate or 
                        read.is_secondary or 
                        read.is_supplementary):
                        total_skipped += 1
                        continue
                    
                    # For paired-end, only process read1
                    if read.is_paired and not read.is_read1:
                        continue
                    
                    # Create fragment
                    fragment = self._create_fragment_from_read(read)
                    if not fragment:
                        total_skipped += 1
                        continue
                    
                    # Process fragment
                    self._process_fragment(fragment)
                    chunk_fragments.append(fragment)
                    total_processed += 1
                    
                    # Process chunk when full
                    if len(chunk_fragments) >= chunk_size:
                        self._write_chunk_to_csv(chunk_fragments, fragment_csv, csv_header_written)
                        if not csv_header_written:
                            csv_header_written = True
                        
                        chunk_count += 1
                        logger.info(f"Processed chunk {chunk_count}: {len(chunk_fragments)} fragments")
                        chunk_fragments = []  # Clear memory
                        
                        # Update progress
                        progress_bar.set_postfix({
                            'chunks': chunk_count,
                            'total_processed': total_processed,
                            'skipped': total_skipped
                        })
                
                # Process final chunk
                if chunk_fragments:
                    self._write_chunk_to_csv(chunk_fragments, fragment_csv, csv_header_written)
                    chunk_count += 1
                    logger.info(f"Processed final chunk {chunk_count}: {len(chunk_fragments)} fragments")
                
                progress_bar.close()
        
        except Exception as e:
            logger.error(f"Error processing BAM file in chunks: {e}")
            raise
        
        logger.info(f"Chunked processing complete: {total_processed:,} total fragments in {chunk_count} chunks")
        
        return {
            'total_fragments': total_processed,
            'total_skipped': total_skipped,
            'chunks_processed': chunk_count,
            'output_file': str(fragment_csv)
        }
    
    def _write_chunk_to_csv(self, fragments: List[OptimizedSalivaFragment], 
                           csv_file: Path, header_written: bool):
        """Write comprehensive fragment chunk to CSV file"""
        mode = 'a' if header_written else 'w'
        
        with open(csv_file, mode, newline='') as f:
            writer = csv.writer(f)
            
            # Write comprehensive header only once
            if not header_written:
                writer.writerow([
                    'read_id', 'chromosome', 'start', 'end', 'length', 'strand',
                    'mapping_quality', 'gc_content', 
                    'predicted_nuclease', 'predicted_5p_nuclease', 'predicted_3p_nuclease', 'nuclease_confidence',
                    'has_g4', 'has_z_dna', 'has_hairpin', 'has_i_motif', 'has_triplex',
                    'mean_mgw', 'mean_helix_twist', 'mean_roll', 'mean_propeller_twist', 'mean_slide', 'mean_shift',
                    'bendability', 'flexibility_score', 'thermal_stability', 'nuclease_accessibility',
                    'damage_score', 'oxidation_count', 'deamination_count', 'abasic_count',
                    'in_oncogenic_region', 'oncogene', 
                    'saliva_nuclease_enriched', 'potential_bacterial_origin'
                ])
            
            # Write comprehensive data
            for frag in fragments:
                writer.writerow([
                    frag.read_id, frag.chromosome, frag.start, frag.end, frag.length, frag.strand,
                    frag.mapping_quality, frag.gc_content,
                    frag.predicted_nuclease, frag.predicted_5p_nuclease, frag.predicted_3p_nuclease, frag.nuclease_confidence,
                    frag.has_g4, frag.has_z_dna, frag.has_hairpin, frag.has_i_motif, frag.has_triplex,
                    frag.mean_mgw, frag.mean_helix_twist, frag.mean_roll, frag.mean_propeller_twist, frag.mean_slide, frag.mean_shift,
                    frag.bendability, frag.flexibility_score, frag.thermal_stability, frag.nuclease_accessibility,
                    frag.damage_score, frag.oxidation_count, frag.deamination_count, frag.abasic_count,
                    frag.in_oncogenic_region, frag.oncogene,
                    frag.saliva_nuclease_enriched, frag.potential_bacterial_origin
                ])

    def process_bam_file(self, bam_path: str, max_fragments: Optional[int] = None) -> List[OptimizedSalivaFragment]:
        """Process BAM file with optimizations"""
        logger.info(f"Processing BAM file: {bam_path}")
        start_time = time.time()
        
        fragments = []
        processed_count = 0
        skipped_count = 0
        
        try:
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                # Get total reads for progress bar (if possible)
                try:
                    total_reads = bam.count()
                    logger.info(f"Total reads in BAM: {total_reads:,}")
                except:
                    total_reads = None
                
                # Process reads with progress bar
                progress_bar = tqdm(
                    bam.fetch(), 
                    total=total_reads,
                    desc="Processing reads",
                    unit="reads"
                )
                
                for read in progress_bar:
                    # Quick filters
                    if (read.is_unmapped or 
                        read.mapping_quality < CONFIG['min_mapping_quality'] or
                        read.is_duplicate or 
                        read.is_secondary or 
                        read.is_supplementary):
                        skipped_count += 1
                        continue
                    
                    # For paired-end, only process read1
                    if read.is_paired and not read.is_read1:
                        continue
                    
                    # Create fragment
                    fragment = self._create_fragment_from_read(read)
                    if not fragment:
                        skipped_count += 1
                        continue
                    
                    # Process fragment
                    self._process_fragment(fragment)
                    
                    fragments.append(fragment)
                    processed_count += 1
                    
                    # Update progress
                    if processed_count % CONFIG['progress_update_interval'] == 0:
                        progress_bar.set_postfix({
                            'processed': processed_count,
                            'skipped': skipped_count
                        })
                    
                    # Check limits
                    if max_fragments and processed_count >= max_fragments:
                        logger.info(f"Reached fragment limit: {max_fragments}")
                        break
                    
                    # Memory management - only stop if limit is set
                    if (CONFIG['max_fragments_in_memory'] is not None and 
                        processed_count >= CONFIG['max_fragments_in_memory']):
                        logger.warning(f"Reached memory limit, stopping at {processed_count} fragments")
                        break
                
                progress_bar.close()
        
        except Exception as e:
            logger.error(f"Error processing BAM file: {e}")
            raise
        
        elapsed_time = time.time() - start_time
        logger.info(f"Processing complete: {processed_count:,} fragments in {elapsed_time:.1f}s")
        logger.info(f"Processing rate: {processed_count/elapsed_time:.1f} fragments/second")
        logger.info(f"Skipped reads: {skipped_count:,}")
        
        return fragments
    
    def _create_fragment_from_read(self, read) -> Optional[OptimizedSalivaFragment]:
        """Create fragment with minimal overhead"""
        # Get fragment coordinates
        if read.is_paired and read.is_proper_pair and not read.mate_is_unmapped:
            frag_start = min(read.reference_start, read.next_reference_start)
            frag_length = abs(read.template_length)
        else:
            frag_start = read.reference_start
            frag_length = read.reference_end - read.reference_start
        
        # Size filter
        if not (CONFIG['min_fragment_length'] <= frag_length <= CONFIG['max_fragment_length']):
            return None
        
        frag_end = frag_start + frag_length
        
        # Extract sequence (only if reference available)
        sequence = None
        gc_content = None
        if self.reference:
            try:
                sequence = self.reference.fetch(
                    read.reference_name, frag_start, frag_end
                ).upper()
                if sequence:
                    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
            except:
                pass  # Skip if sequence extraction fails
        
        # Create fragment
        fragment = OptimizedSalivaFragment(
            read_id=read.query_name,
            chromosome=read.reference_name,
            start=frag_start,
            end=frag_end,
            length=frag_length,
            strand='-' if read.is_reverse else '+',
            mapping_quality=read.mapping_quality,
            sequence=sequence,
            gc_content=gc_content
        )
        
        # Quick oncogenic region check
        if fragment.chromosome in self.oncogenic_lookup:
            for start, end, gene in self.oncogenic_lookup[fragment.chromosome]:
                if start <= fragment.start <= end or start <= fragment.end <= end:
                    fragment.in_oncogenic_region = True
                    fragment.oncogene = gene
                    break
        
        return fragment
    
    def _process_fragment(self, fragment: OptimizedSalivaFragment):
        """Process fragment with enabled analyzers only"""
        if not fragment.sequence:
            return
        
        # Structure detection
        if self.structure_detector and CONFIG['enable_structure_detection']:
            structures = self.structure_detector.detect_structures(fragment.sequence)
            fragment.has_g4 = structures['has_g4']
            fragment.has_z_dna = structures['has_z_dna']
            fragment.has_hairpin = structures['has_hairpin']
        
        # Damage analysis
        if self.damage_analyzer and CONFIG['enable_damage_analysis']:
            damage_results = self.damage_analyzer.analyze_damage(fragment.sequence)
            fragment.damage_score = damage_results['damage_score']
            fragment.oxidation_count = damage_results['oxidation_count']
            fragment.deamination_count = damage_results['deamination_count']
        
        # Nuclease prediction
        nuclease_5p, conf_5p, nuclease_3p, conf_3p = self.nuclease_predictor.predict_nuclease(fragment.sequence)
        fragment.predicted_5p_nuclease = nuclease_5p
        fragment.predicted_3p_nuclease = nuclease_3p
        # Use the higher confidence prediction as the main prediction
        if conf_5p >= conf_3p:
            fragment.predicted_nuclease = nuclease_5p
            fragment.nuclease_confidence = conf_5p
        else:
            fragment.predicted_nuclease = nuclease_3p
            fragment.nuclease_confidence = conf_3p    


            # Check if saliva-enriched
            if fragment.predicted_nuclease and COMPREHENSIVE_NUCLEASE_SIGNATURES.get(fragment.predicted_nuclease, {}).get('saliva_enriched', False):
                fragment.saliva_nuclease_enriched = True

# ============================================================================
# OPTIMIZED ANALYSIS COORDINATOR
# ============================================================================

class OptimizedAnalysisCoordinator:
    """Fast analysis coordinator"""
    
    def __init__(self):
        pass
    
    def run_analysis(self, fragments: List[OptimizedSalivaFragment]) -> Dict[str, Any]:
        """Run comprehensive analysis"""
        logger.info("Running comprehensive analysis...")
        start_time = time.time()
        
        results = {
            'fragment_statistics': self._calculate_fragment_stats(fragments),
            'nuclease_analysis': self._analyze_nucleases(fragments),
            'structure_analysis': self._analyze_structures(fragments),
            'damage_analysis': self._analyze_damage(fragments),
            'shape_analysis': self._analyze_shape_features(fragments),
            'saliva_features': self._analyze_saliva_features(fragments),
            'summary': {}
        }
        
        # Generate summary
        results['summary'] = self._generate_summary(results, fragments)
        
        elapsed_time = time.time() - start_time
        logger.info(f"Comprehensive analysis complete in {elapsed_time:.1f}s")
        
        return results
    
    def _calculate_fragment_stats(self, fragments: List[OptimizedSalivaFragment]) -> Dict:
        """Calculate basic statistics efficiently"""
        if not fragments:
            return {}
        
        lengths = np.array([f.length for f in fragments])
        gc_contents = np.array([f.gc_content for f in fragments if f.gc_content is not None])
        
        return {
            'total_fragments': len(fragments),
            'length_stats': {
                'mean': float(np.mean(lengths)),
                'median': float(np.median(lengths)),
                'std': float(np.std(lengths)),
                'min': int(np.min(lengths)),
                'max': int(np.max(lengths))
            },
            'gc_stats': {
                'mean': float(np.mean(gc_contents)) if len(gc_contents) > 0 else 0,
                'std': float(np.std(gc_contents)) if len(gc_contents) > 0 else 0
            } if len(gc_contents) > 0 else {}
        }
    
    def _analyze_nucleases(self, fragments: List[OptimizedSalivaFragment]) -> Dict:
        """Analyze nuclease patterns"""
        nuclease_counts = Counter()
        confidence_scores = defaultdict(list)
        
        for fragment in fragments:
            if fragment.predicted_nuclease:
                nuclease_counts[fragment.predicted_nuclease] += 1
                confidence_scores[fragment.predicted_nuclease].append(fragment.nuclease_confidence)
        
        results = {}
        total_fragments = len(fragments)
        
        for nuclease, count in nuclease_counts.items():
            results[nuclease] = {
                'count': count,
                'percentage': count / total_fragments * 100,
                'mean_confidence': np.mean(confidence_scores[nuclease])
            }
        
        return results
    
    def _analyze_structures(self, fragments: List[OptimizedSalivaFragment]) -> Dict:
        """Analyze comprehensive structural patterns"""
        structure_counts = {
            'g4': sum(1 for f in fragments if f.has_g4),
            'z_dna': sum(1 for f in fragments if f.has_z_dna),
            'hairpin': sum(1 for f in fragments if f.has_hairpin),
            'i_motif': sum(1 for f in fragments if f.has_i_motif),
            'triplex': sum(1 for f in fragments if f.has_triplex)
        }
        
        total_fragments = len(fragments)
        
        # Multiple structure analysis
        multi_structure_count = sum(1 for f in fragments 
                                  if sum([f.has_g4, f.has_z_dna, f.has_hairpin, f.has_i_motif, f.has_triplex]) > 1)
        
        return {
            'counts': structure_counts,
            'percentages': {k: v / total_fragments * 100 for k, v in structure_counts.items()},
            'fragments_with_structures': sum(1 for f in fragments 
                                           if f.has_g4 or f.has_z_dna or f.has_hairpin or f.has_i_motif or f.has_triplex),
            'fragments_with_multiple_structures': multi_structure_count,
            'structure_diversity_index': len([k for k, v in structure_counts.items() if v > 0])
        }
    
    def _analyze_damage(self, fragments: List[OptimizedSalivaFragment]) -> Dict:
        """Analyze comprehensive damage patterns"""
        damage_scores = np.array([f.damage_score for f in fragments])
        oxidation_counts = np.array([f.oxidation_count for f in fragments])
        deamination_counts = np.array([f.deamination_count for f in fragments])
        abasic_counts = np.array([f.abasic_count for f in fragments])
        
        high_damage_threshold = 0.1
        high_damage_count = np.sum(damage_scores > high_damage_threshold)
        
        return {
            'mean_damage_score': float(np.mean(damage_scores)),
            'median_damage_score': float(np.median(damage_scores)),
            'damage_score_std': float(np.std(damage_scores)),
            'high_damage_count': int(high_damage_count),
            'high_damage_percentage': float(high_damage_count / len(fragments) * 100),
            'oxidation_stats': {
                'mean': float(np.mean(oxidation_counts)),
                'max': int(np.max(oxidation_counts)),
                'total': int(np.sum(oxidation_counts))
            },
            'deamination_stats': {
                'mean': float(np.mean(deamination_counts)),
                'max': int(np.max(deamination_counts)),
                'total': int(np.sum(deamination_counts))
            },
            'abasic_stats': {
                'mean': float(np.mean(abasic_counts)),
                'max': int(np.max(abasic_counts)),
                'total': int(np.sum(abasic_counts))
            }
        }
    
    def _analyze_shape_features(self, fragments: List[OptimizedSalivaFragment]) -> Dict:
        """Analyze DNA shape features if available"""
        if not CONFIG['enable_shape_analysis']:
            return {'shape_analysis_enabled': False}
        
        # Extract shape values
        mgw_values = [f.mean_mgw for f in fragments if f.mean_mgw > 0]
        twist_values = [f.mean_helix_twist for f in fragments if f.mean_helix_twist > 0]
        roll_values = [f.mean_roll for f in fragments]
        bendability_values = [f.bendability for f in fragments if f.bendability > 0]
        flexibility_values = [f.flexibility_score for f in fragments if f.flexibility_score > 0]
        
        if not mgw_values:
            return {'shape_analysis_enabled': True, 'insufficient_data': True}
        
        return {
            'shape_analysis_enabled': True,
            'mgw_stats': {
                'mean': np.mean(mgw_values),
                'std': np.std(mgw_values),
                'min': np.min(mgw_values),
                'max': np.max(mgw_values)
            },
            'twist_stats': {
                'mean': np.mean(twist_values) if twist_values else 0,
                'std': np.std(twist_values) if twist_values else 0
            },
            'bendability_stats': {
                'mean': np.mean(bendability_values) if bendability_values else 0,
                'std': np.std(bendability_values) if bendability_values else 0,
                'high_bendability_count': sum(1 for b in bendability_values if b > 0.7)
            },
            'flexibility_stats': {
                'mean': np.mean(flexibility_values) if flexibility_values else 0,
                'std': np.std(flexibility_values) if flexibility_values else 0,
                'high_flexibility_count': sum(1 for f in flexibility_values if f > 0.7)
            }
        }
    
    def _analyze_saliva_features(self, fragments: List[OptimizedSalivaFragment]) -> Dict:
        """Analyze saliva-specific features"""
        saliva_enriched_count = sum(1 for f in fragments if f.saliva_nuclease_enriched)
        oncogenic_count = sum(1 for f in fragments if f.in_oncogenic_region)
        
        # Size distribution
        lengths = [f.length for f in fragments]
        size_categories = {
            'ultra_short': sum(1 for l in lengths if l < 50),
            'short': sum(1 for l in lengths if 50 <= l < 100),
            'medium': sum(1 for l in lengths if 100 <= l < 200),
            'long': sum(1 for l in lengths if 200 <= l < 400),
            'very_long': sum(1 for l in lengths if l >= 400)
        }
        
        return {
            'saliva_nuclease_enriched': {
                'count': saliva_enriched_count,
                'percentage': saliva_enriched_count / len(fragments) * 100
            },
            'oncogenic_regions': {
                'count': oncogenic_count,
                'percentage': oncogenic_count / len(fragments) * 100
            },
            'size_distribution': size_categories
        }
    
    def _generate_summary(self, results: Dict, fragments: List[OptimizedSalivaFragment]) -> Dict:
        """Generate analysis summary"""
        key_findings = []
        
        # Check for high saliva nuclease signature
        if 'saliva_features' in results:
            saliva_pct = results['saliva_features']['saliva_nuclease_enriched']['percentage']
            if saliva_pct > 30:
                key_findings.append(f"High saliva nuclease signature: {saliva_pct:.1f}%")
        
        # Check for damage
        if 'damage_analysis' in results:
            damage_pct = results['damage_analysis']['high_damage_percentage']
            if damage_pct > 20:
                key_findings.append(f"Elevated DNA damage: {damage_pct:.1f}%")
        
        # Top nuclease
        if 'nuclease_analysis' in results and results['nuclease_analysis']:
            top_nuclease = max(results['nuclease_analysis'].items(), 
                             key=lambda x: x[1]['count'])
            key_findings.append(f"Top nuclease: {top_nuclease[0]} ({top_nuclease[1]['percentage']:.1f}%)")
        
        return {
            'analysis_timestamp': datetime.now().isoformat(),
            'total_fragments': len(fragments),
            'key_findings': key_findings,
            'processing_config': {
                'structure_detection': CONFIG['enable_structure_detection'],
                'damage_analysis': CONFIG['enable_damage_analysis'],
                'nuclease_prediction': CONFIG['enable_nuclease_prediction'],
                'shape_analysis': CONFIG['enable_shape_analysis']
            }
        }
    
    def export_results_csv(self, fragments: List[OptimizedSalivaFragment], 
                          results: Dict, output_dir: Path):
        """Export comprehensive results to CSV efficiently"""
        logger.info("Exporting comprehensive results to CSV...")
        
        # Main fragment data with all comprehensive fields
        fragment_file = output_dir / 'comprehensive_saliva_fragments.csv'
        with open(fragment_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Comprehensive header
            writer.writerow([
                'read_id', 'chromosome', 'start', 'end', 'length', 'strand',
                'mapping_quality', 'gc_content', 
                'predicted_nuclease', 'predicted_5p_nuclease', 'predicted_3p_nuclease', 'nuclease_confidence',
                'has_g4', 'has_z_dna', 'has_hairpin', 'has_i_motif', 'has_triplex',
                'mean_mgw', 'mean_helix_twist', 'mean_roll', 'mean_propeller_twist', 'mean_slide', 'mean_shift',
                'bendability', 'flexibility_score', 'thermal_stability', 'nuclease_accessibility',
                'damage_score', 'oxidation_count', 'deamination_count', 'abasic_count',
                'in_oncogenic_region', 'oncogene', 
                'saliva_nuclease_enriched', 'potential_bacterial_origin'
            ])
            
            # Comprehensive data
            for frag in fragments:
                writer.writerow([
                    frag.read_id, frag.chromosome, frag.start, frag.end, frag.length, frag.strand,
                    frag.mapping_quality, frag.gc_content,
                    frag.predicted_nuclease, frag.predicted_5p_nuclease, frag.predicted_3p_nuclease, frag.nuclease_confidence,
                    frag.has_g4, frag.has_z_dna, frag.has_hairpin, frag.has_i_motif, frag.has_triplex,
                    frag.mean_mgw, frag.mean_helix_twist, frag.mean_roll, frag.mean_propeller_twist, frag.mean_slide, frag.mean_shift,
                    frag.bendability, frag.flexibility_score, frag.thermal_stability, frag.nuclease_accessibility,
                    frag.damage_score, frag.oxidation_count, frag.deamination_count, frag.abasic_count,
                    frag.in_oncogenic_region, frag.oncogene,
                    frag.saliva_nuclease_enriched, frag.potential_bacterial_origin
                ])
        
        # Nuclease analysis summary
        nuclease_file = output_dir / 'nuclease_analysis_comprehensive.csv'
        with open(nuclease_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['nuclease', 'fragment_count', 'percentage', 'mean_confidence', 'saliva_enriched'])
            
            nuclease_counts = defaultdict(int)
            nuclease_confidences = defaultdict(list)
            
            for frag in fragments:
                if frag.predicted_nuclease:
                    nuclease_counts[frag.predicted_nuclease] += 1
                    nuclease_confidences[frag.predicted_nuclease].append(frag.nuclease_confidence)
            
            total_fragments = len(fragments)
            for nuclease, count in nuclease_counts.items():
                is_saliva_enriched = COMPREHENSIVE_NUCLEASE_SIGNATURES.get(nuclease, {}).get('saliva_enriched', False)
                mean_conf = np.mean(nuclease_confidences[nuclease]) if nuclease_confidences[nuclease] else 0
                writer.writerow([
                    nuclease, count, count / total_fragments * 100, mean_conf, is_saliva_enriched
                ])
        
        # Shape analysis summary
        if CONFIG['enable_shape_analysis']:
            shape_file = output_dir / 'shape_analysis_summary.csv'
            with open(shape_file, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['metric', 'mean', 'std', 'min', 'max', 'q25', 'q50', 'q75'])
                
                shape_metrics = {
                    'mgw': [f.mean_mgw for f in fragments if f.mean_mgw > 0],
                    'helix_twist': [f.mean_helix_twist for f in fragments if f.mean_helix_twist > 0],
                    'roll': [f.mean_roll for f in fragments],
                    'propeller_twist': [f.mean_propeller_twist for f in fragments],
                    'bendability': [f.bendability for f in fragments],
                    'flexibility': [f.flexibility_score for f in fragments],
                    'thermal_stability': [f.thermal_stability for f in fragments],
                    'nuclease_accessibility': [f.nuclease_accessibility for f in fragments]
                }
                
                for metric, values in shape_metrics.items():
                    if values:
                        writer.writerow([
                            metric, np.mean(values), np.std(values), np.min(values), np.max(values),
                            np.percentile(values, 25), np.percentile(values, 50), np.percentile(values, 75)
                        ])
        
        # Summary results
        summary_file = output_dir / 'comprehensive_analysis_summary.csv'
        with open(summary_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['metric', 'value'])
            
            # Flatten results for CSV
            def flatten_dict(d, prefix=''):
                items = []
                for k, v in d.items():
                    if isinstance(v, dict):
                        items.extend(flatten_dict(v, f"{prefix}{k}_"))
                    else:
                        items.append((f"{prefix}{k}", v))
                return items
            
            for metric, value in flatten_dict(results):
                writer.writerow([metric, value])

# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    """Comprehensive main function with full analysis options"""
    parser = argparse.ArgumentParser(
        description='Comprehensive Saliva cfDNA Analysis Pipeline - Full feature analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full comprehensive analysis (all features enabled)
  python script.py -i input.bam -r reference.fasta
  
  # Large file with chunked processing
  python script.py -i large_file.bam -r reference.fasta --chunked-processing --chunk-size 100000
  
  # Quick analysis (disable expensive features)
  python script.py -i input.bam --disable-shape --chunk-size 20000
  
  # Test run with limited fragments
  python script.py -i input.bam -r reference.fasta --max-fragments 10000
        """
    )
    parser.add_argument('-i', '--input', required=True, help='Input BAM file')
    parser.add_argument('-r', '--reference', help='Reference genome FASTA (optional but recommended)')
    parser.add_argument('-o', '--output', default='comprehensive_saliva_analysis', help='Output directory')
    parser.add_argument('--max-fragments', type=int, help='Maximum fragments to analyze (for testing)')
    parser.add_argument('--chunk-size', type=int, default=50000, help='Processing chunk size')
    parser.add_argument('--chunked-processing', action='store_true', 
                       help='Use chunked processing for very large files (saves memory)')
    
    # Analysis control options
    parser.add_argument('--disable-structures', action='store_true', help='Disable structure detection')
    parser.add_argument('--disable-damage', action='store_true', help='Disable damage analysis')
    parser.add_argument('--disable-nuclease', action='store_true', help='Disable nuclease prediction')
    parser.add_argument('--disable-shape', action='store_true', help='Disable DNA shape analysis (faster)')
    
    # Performance options
    parser.add_argument('--fast-mode', action='store_true', 
                       help='Enable fast mode (disables shape analysis, larger chunks)')
    parser.add_argument('--workers', type=int, help='Number of worker processes')
    
    args = parser.parse_args()
    
    # Apply fast mode settings
    if args.fast_mode:
        CONFIG['enable_shape_analysis'] = False
        CONFIG['chunk_size'] = 100000
        logger.info("Fast mode enabled: shape analysis disabled, larger chunks")
    
    # Update config based on arguments
    if args.chunk_size:
        CONFIG['chunk_size'] = args.chunk_size
    if args.workers:
        CONFIG['max_workers'] = min(args.workers, mp.cpu_count())
    if args.disable_structures:
        CONFIG['enable_structure_detection'] = False
    if args.disable_damage:
        CONFIG['enable_damage_analysis'] = False
    if args.disable_nuclease:
        CONFIG['enable_nuclease_prediction'] = False
    if args.disable_shape:
        CONFIG['enable_shape_analysis'] = False
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("Starting comprehensive saliva cfDNA analysis")
    logger.info(f"Analysis configuration:")
    logger.info(f"  Shape analysis: {CONFIG['enable_shape_analysis']}")
    logger.info(f"  Structure detection: {CONFIG['enable_structure_detection']}")
    logger.info(f"  Damage analysis: {CONFIG['enable_damage_analysis']}")
    logger.info(f"  Nuclease prediction: {CONFIG['enable_nuclease_prediction']}")
    logger.info(f"  Chunk size: {CONFIG['chunk_size']:,}")
    logger.info(f"  Chunked processing: {args.chunked_processing}")
    
    total_start_time = time.time()
    
    try:
        # Initialize processor
        processor = OptimizedFragmentProcessor(args.reference)
        
        # Choose processing method
        if args.chunked_processing:
            # Chunked processing for large files
            logger.info("Using chunked processing mode for large files")
            results = processor.process_bam_file_chunked(args.input, output_dir, args.chunk_size)
            
            # For chunked processing, we don't do full analysis
            print("\n" + "="*80)
            print("CHUNKED PROCESSING COMPLETE")
            print("="*80)
            print(f"Total fragments processed: {results['total_fragments']:,}")
            print(f"Chunks processed: {results['chunks_processed']}")
            print(f"Results saved to: {results['output_file']}")
            
        else:
            # Standard processing (loads all fragments into memory)
            fragments = processor.process_bam_file(args.input, args.max_fragments)
            
            if not fragments:
                logger.error("No fragments were processed!")
                return 1
            
            # Run analysis
            coordinator = OptimizedAnalysisCoordinator()
            results = coordinator.run_analysis(fragments)
            
            # Save JSON results
            json_file = output_dir / 'comprehensive_analysis_results.json'
            with open(json_file, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            
            # Export CSV
            coordinator.export_results_csv(fragments, results, output_dir)
            
            # Print results
            total_time = time.time() - total_start_time
            
            print("\n" + "="*80)
            print("COMPREHENSIVE SALIVA cfDNA ANALYSIS COMPLETE")
            print("="*80)
            print(f"Total runtime: {total_time:.1f} seconds")
            print(f"Fragments analyzed: {len(fragments):,}")
            print(f"Processing rate: {len(fragments)/total_time:.1f} fragments/second")
            
            # Analysis summary
            enabled_analyses = []
            if CONFIG['enable_shape_analysis']:
                enabled_analyses.append("DNA Shape")
            if CONFIG['enable_structure_detection']:
                enabled_analyses.append("Structure Detection")
            if CONFIG['enable_damage_analysis']:
                enabled_analyses.append("Damage Analysis")
            if CONFIG['enable_nuclease_prediction']:
                enabled_analyses.append("Nuclease Prediction")
            
            print(f"Enabled analyses: {', '.join(enabled_analyses)}")
            
            # Key statistics
            if 'fragment_statistics' in results:
                stats = results['fragment_statistics']
                if 'length_stats' in stats:
                    print(f"Mean fragment length: {stats['length_stats']['mean']:.1f} bp")
                    print(f"Fragment range: {stats['length_stats']['min']}-{stats['length_stats']['max']} bp")
            
            # Shape analysis summary
            if 'shape_analysis' in results and results['shape_analysis'].get('shape_analysis_enabled'):
                if 'mgw_stats' in results['shape_analysis']:
                    mgw_stats = results['shape_analysis']['mgw_stats']
                    print(f"Mean minor groove width: {mgw_stats['mean']:.2f} ± {mgw_stats['std']:.2f} Å")
                if 'bendability_stats' in results['shape_analysis']:
                    bend_stats = results['shape_analysis']['bendability_stats']
                    print(f"High bendability fragments: {bend_stats.get('high_bendability_count', 0):,}")
            
            # Structure analysis summary
            if 'structure_analysis' in results:
                struct_counts = results['structure_analysis']['counts']
                total_structures = sum(struct_counts.values())
                if total_structures > 0:
                    print(f"Total structural motifs detected: {total_structures:,}")
                    top_structure = max(struct_counts.items(), key=lambda x: x[1])
                    print(f"Most common structure: {top_structure[0]} ({top_structure[1]:,} fragments)")
            
            # Nuclease analysis summary
            if 'nuclease_analysis' in results and results['nuclease_analysis']:
                top_nuclease = max(results['nuclease_analysis'].items(), 
                                 key=lambda x: x[1]['count'])
                print(f"Top nuclease: {top_nuclease[0]} ({top_nuclease[1]['count']:,} fragments, {top_nuclease[1]['percentage']:.1f}%)")
            
            # Key findings
            if 'summary' in results and 'key_findings' in results['summary']:
                print("\nKey Findings:")
                for finding in results['summary']['key_findings']:
                    print(f"  • {finding}")
            
            print(f"\nResults saved to: {output_dir}")
            print("Output files:")
            # print(f"  • comprehensive_saliva_fragments.csv - Main fragment data")
            print(f"  • nuclease_analysis_comprehensive.csv - Nuclease signatures")
            if CONFIG['enable_shape_analysis']:
                print(f"  • shape_analysis_summary.csv - DNA shape statistics")
            print(f"  • comprehensive_analysis_summary.csv - Overall summary")
            print(f"  • comprehensive_analysis_results.json - Full results")
        
        print("="*80)
        return 0
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())