"""
Author: Andrew Harris
Python version 3.8
"""
import logging
import os
from pathlib import Path
import random

import pandas as pd
import vcfpy

# ---------------- Helper Functions ----------------
def roundup(x, window_size):
    """
    This function takes in a position and rounds it up to the nearest 100,000th position.
    @param x: Position of record
    @return:Position rounded up to the nearest 100,000th.
    """
    return x if x % window_size == 0 else x + window_size - x % window_size


def collect_introgression_info(sample_name, chromosome, position, output_dir, project, window_size_str):
    """This function takes in the sample name of interest, the chromosome, and the position of the record rounded to
    the nearest window size (i.e. 100-kb).
    return: True if introgressed"""
    file_pathway = f'{output_dir}/{project}_{window_size_str}/zscore_boolean/{sample_name}.zscore_boolean.tsv'
    df = pd.read_csv(file_pathway, sep='\t', index_col=[0])
    is_introgressed = df.at[position, str(chromosome)]
    if is_introgressed:
        return True
    else:
        return False


def get_samples(excel_file):
    """Collect cow and bison sample names and return them in a list"""
    cow_samples = []
    bison_samples = []

    sample_info_file = pd.read_excel(
        excel_file,
        usecols=["Sample", "Subspecies_Reference_Number"],
    )

    for _, row in sample_info_file.iterrows():
        sample_name = list(row)[0]
        species_num = list(row)[1]
        if species_num == 0:
            cow_samples.append(sample_name)
        elif species_num == 1:
            bison_samples.append(sample_name)
        else:
            # Outgroup number
            continue

    return cow_samples, bison_samples

# ----------------- Main Function -----------------
def alter_original_vcf(
    input_directory,
    output_directory,
    project,
    window_size_str,
    window_size_int,
    original_file,
    chromosome_name,
    sample_reference_info,
    logger,
):
    """
    This function will take in the introgression_info_dict that was returned by collect_introgression_info
    and will alter the original vcf file to reflect the introgression events we have identified. All samples for
    a given record that do not meet the z-score threshold will be set homozygous for the alternative
    """

    wp = f'{output_directory}/{project}_{window_size_str}/updated_vcf_output/'
    writer_outdir = Path(wp)
    writer_outdir.mkdir(parents=True, exist_ok=True)

    out_file_name = writer_outdir / f'{chromosome_name}_altered_haplotypes.vcf'
    reader = vcfpy.Reader.from_path(original_file)
    reader.header.add_info_line(
        vcfpy.OrderedDict(
            [
                ('ID', 'INTROGRESSED'),
                ("Number", 1),
                ("Type", "String"),
                ('Description', 'Indication if record is introgressed or not')
            ]
        )
    )
    writer = vcfpy.Writer.from_path(
        out_file_name,
        reader.header,
    )

    # Tracking variables
    current_rec = 1
    total_records_altered = 0

    # Collect sample names by species
    cow_samples, bison_samples = get_samples(sample_reference_info)

    # mono+tri-allelic sites
    tri_alelle_df = pd.DataFrame(columns=["Chromsome", "Position", "Sample", "Genotype", "Introgressed"])
    single_haplotype_df = pd.DataFrame(columns=["Chromsome", "Position", "Sample", "Genotype", "Introgressed"])

    tri_idx = 0
    sh_idx = 0

    # Iterate through all records
    for record in reader:
        # Set Introgression to False. Only changed to True if sample in record is introgressed. Does not
        # indicate which sample though.
        record.INFO['INTROGRESSED'] = 'F'

        init_record = str(record)
        current_rec += 1

        # Record column values
        chromosome = record.CHROM
        rounded_position = roundup(record.POS, window_size_int)

        # If the record is not a snv, ignore and continue
        if not record.is_snv():
            continue

        for call in record.calls:
            # Ignores . calls
            # We are also only interest in het sites
            if call.is_het:
                # Call breakdown for easy reference
                sample = call.sample
                genotype = call.data['GT']
                allele_1 = int(genotype[0])  # Type=int
                allele_2 = int(genotype[2])  # Type=int
                # Check if sample is bison, if so call collect_introgression_info(sample, CHROM, rounded POS)
                if sample in bison_samples:
                    # Check to see if the sample is introgressed
                    is_introgressed = collect_introgression_info(
                        sample,
                        chromosome,
                        rounded_position,
                        output_directory,
                        project,
                        window_size_str,
                    )
                    if is_introgressed:
                        record.INFO['INTROGRESSED'] = f"T-{sample}"
                        # Check for introgressed Triplet
                        # Randomly choose one of the alternatives to represent
                        if allele_1 != 0:
                            if allele_2 != 0:
                                # check to see if triplet is found in cattle
                                for cow_call in record.calls:
                                    if cow_call.sample in cow_samples:
                                        if cow_call.data['GT'] == '.':
                                            continue
                                        elif '0' in cow_call.data['GT']:
                                            # Skipped: cow genotype contained reference allele
                                            continue
                                        elif cow_call.data['GT'] == call.data["GT"]:
                                            # Skipped because cow shows triplet genotype
                                            continue
                                        else:
                                            # There are some genotypes that only have 1 int not 0/1
                                            if type(call.data['GT']) == int:
                                                single_haplotype_df.at[sh_idx, 'Chromosome'].append(record.CHROM)
                                                single_haplotype_df.at[sh_idx, 'Position'].append(record.POS)
                                                single_haplotype_df.at[sh_idx, 'Sample'].append(sample)
                                                single_haplotype_df.at[sh_idx, 'Genotype'].append(genotype)
                                                single_haplotype_df.at[sh_idx, 'Introgressed'].append('True')
                                                sh_idx += 1
                                                continue
                                            else:
                                                tri_alelle_df.at[tri_idx, 'Chromosome'].append(record.CHROM)
                                                tri_alelle_df.at[tri_idx, 'Position'].append(record.POS)
                                                tri_alelle_df.at[tri_idx, 'Sample'].append(sample)
                                                tri_alelle_df.at[tri_idx, 'Genotype'].append(genotype)
                                                tri_alelle_df.at[tri_idx, 'Introgressed'].append('True')
                                                tri_idx += 1
                                                # Randomly choose one of the alleles
                                                call.data['GT'] = random.choice([allele_1, allele_2])
                                                continue
                                    else:
                                        continue
                                continue
                            else:
                                call.data['GT'] = '0/0'
                                logger.debug(call)
                                pass
                        # Reads that end up here are not triplets, but heterozygotes (ref/alt)
                        else:
                            call.data['GT'] = call.data['GT'] = '0/0'
                            logger.debug(call)
                            continue
                    # Identify non-introgressed heterozygotes and make them homozygous alternative
                    # If allele 1 is reference (0), assign genotype to homozygous alelle 2.
                    elif allele_1 == 0:
                        if record.INFO['INTROGRESSED'] == 'T':
                            call.data['GT'] = f"{allele_2}/{allele_2}"
                        else:
                            record.INFO['INTROGRESSED'] = 'F'
                            call.data['GT'] = f"{allele_2}/{allele_2}"
                        continue
                    # If allele 2 is reference (0), assign genotype to homozygous allele 1.
                    elif allele_2 == 0:
                        if record.INFO['INTROGRESSED'] == 'T':
                            call.data['GT'] = f"{allele_1}/{allele_1}"
                        else:
                            record.INFO['INTROGRESSED'] = 'F'
                            call.data['GT'] = f"{allele_1}/{allele_1}"
                        continue
                    # If both alleles are not 0
                    # Check if it is homozygous alternative
                    elif allele_1 == allele_2:
                        # Do nothing to genotype
                        call.data['INTRO'] = 'no'
                        logger.debug(call.data)
                        logger.debug(call.data['GT'])
                        logger.debug("Homozygous Alternative")
                        continue
                    # Check for non-introgressed Triplets
                    elif allele_1 != 0:
                        if allele_2 != 0:
                            # check to see if triplet is found in cattle
                            for cow_call in record.calls:
                                if cow_call.sample in cow_samples:
                                    # Check is sample is a no call
                                    if cow_call.data['GT'] == '.':
                                        # Skipped because it is a no call
                                        continue
                                    # If cow is contained reference allele
                                    elif '0' in cow_call.data['GT']:
                                        # Skipped cow genotype contained reference allele
                                        continue
                                    elif cow_call.data['GT'] == call.data["GT"]:
                                        # Skipped because cow shows triplet genotype
                                        continue
                                    else:
                                        # If genotype is of type int, this means it only has one allele, which is an issue...
                                        if type(call.data['GT']) == int:
                                            logger.debug("Non-introgressed Single haplotype ISSUE!!!!!!")
                                            logger.debug(f"{record.CHROM}/{record.POS}")
                                            logger.debug(call.data['GT'])
                                            logger.debug(f"{allele_1}/{allele_2}")
                                            single_haplotype_df.at[sh_idx, 'Chromosome'].append(record.CHROM)
                                            single_haplotype_df.at[sh_idx, 'Position'].append(record.POS)
                                            single_haplotype_df.at[sh_idx, 'Sample'].append(sample)
                                            single_haplotype_df.at[sh_idx, 'Genotype'].append(genotype)
                                            single_haplotype_df.at[sh_idx, 'Introgressed'].append('False')
                                            sh_idx += 1
                                            continue
                                        else:  # Found Triplet, record.
                                            tri_alelle_df.at[tri_idx, 'Chromosome'].append(record.CHROM)
                                            tri_alelle_df.at[tri_idx, 'Position'].append(record.POS)
                                            tri_alelle_df.at[tri_idx, 'Sample'].append(sample)
                                            tri_alelle_df.at[tri_idx, 'Genotype'].append(genotype)
                                            tri_alelle_df.at[tri_idx, 'Introgressed'].append('False')
                                            tri_idx += 1
                                            call.data['GT'] = random.choice([allele_1, allele_2])
                                            call.data['INTRO'] = 'no'
                                            continue
                                else:
                                    continue
                            continue
                        else:
                            pass
                    # Otherwise switch genotype to homozygous alternative (bison-like)
                    else:
                        record.INFO['INTROGRESSED'] = 'F'
                        call.data['GT'] = f'{allele_2}/{allele_2}'
                        logger.debug(call.data['GT'])
                        continue
                else:
                    continue
            else:
                continue
        if init_record != record:
            total_records_altered += 1
            writer.write_record(record)
            continue
        else:
            writer.write_record(record)
            continue
    # tri_allele_output_directory
    tra_outpu_dir = Path(f'{output_directory}/{project}_{window_size_str}/updated_vcf_output/tri_allele/')
    # single_haplotype_output_directory
    sh_output_dir = Path(f'{output_directory}/{project}_{window_size_str}/updated_vcf_output/single_haplotype/')
    # Make output dirs if not already created
    tra_outpu_dir.mkdir(parents=True, exist_ok=True)
    sh_output_dir.mkdir(parents=True, exist_ok=True)

    # tri_alelle_df = pd.DataFrame().from_dict(tri_allele_collect)
    tri_allele_output_file = tra_outpu_dir / f"{chromosome_name}_tri_allele_positions.tsv"
    tri_alelle_df.to_csv(tri_allele_output_file, sep='\t')

    # single_haplotype_issue_df = pd.DataFrame().from_dict(single_haplotype_issue_tracking)
    single_haplotype_output_file = sh_output_dir / f"{chromosome_name}_haploid_records.tsv"
    single_haplotype_df.to_csv(single_haplotype_output_file, sep='\t')
    return
