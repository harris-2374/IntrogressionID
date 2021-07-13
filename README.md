# IntrogressionID
IntrogressionID is a pipeline of scripts that takes a VCF file containing variants from a population of introgressed individuals  


# Pipeline:
    1. Allele frequency filter
    2. Windowed SNP counting
    3. Z-score distribution
    4. Plot per-sample heatmaps and manhattan plots
    5. Update original VCF files with introgression data

# Steps:
You can run IntrogressionID step-by-step using one of the six arguments below. Providing more than one argument will
 run all steps in the order they come in the pipeline. Providing no argument will run the entire pipeline start to
  finish. 
 
    (--allele) -- [run allele frequency filtration step]
    (--vcf_check) -- [run secondary variant cross-check step]
    (--count) -- [run windowed SNP counting step]
    (--zscore) -- [run z-score distribution step]
    (--plot) -- [plot per-sample heatmaps and manhattan plots]
    (--alter) -- [update original vcf files]
    
 # Input Arguments:
 
    Required arguments for all steps:
        (--project) (-p) -- [provide a project name to organize data]
        (--window_size) (-w) -- [sliding window size (i.e. 5bp, 5kb, 5mb)]
        
    Step dependednt arguments: (refer to example commands below)
        (--input) (-i)  -- [pathway to directory containing input vcf files broken down per-chromosome]
        (--output) (-o) -- [output directory pathway]
        (--reference) (-r) -- [filepath to excel sheet containing sample species information]
        (--bed) (-b) -- [bed file with chromosome lengths]
        (--threshold) (-t) -- [z-score threshold]
 
# Options:
    (--auto_graph) -- [opens all graphs in the default broswer; not reccomended when running "--heatmaps"]


# Installation and Setup:
    
    
# Usage:

    Example command: python introgressionID.py  [ step ] [ inputs ] [ options ]
    
    You can run the entire pipeline by providing no step argument:
        python introgressionID.py \
        -p Project_Name \
        -w 100kb \
        -i chromosome_data/ \
        -r sample_reference_data.xlsx \
        -b Chromosome_Lengths.bed \
        -t 5.0
        
    
    You can also run each step individually:
    
        python introgressionID.py --allele \
        -i per_chromosome_vcf_files/ \
        -p Project_Name \
        -w 100kb \
        -r sample_reference_data.xlsx
        
        python introgressionID.py --count \
        -p Project_Name \
        -w 100kb \
        -b Chromosome_Lengths.bed \
        
        python introgressionID.py --zscore \
        -p Project_Name \
        -w 100kb \
        -b Chromosome_Lengths.bed \
        -t 5.0
        
        python introgressionID.py --plot \
        -p Project_Name \
        -w 100kb \
        -t 5.0 \
        
        python introgressionID.py --alter \
        -i per_chromosome_vcf_files/ \
        -p Project_Name \
        -w 100kb \
        -r sample_reference_data.xlsx
        
    
