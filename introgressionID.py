import argparse
from datetime import datetime
import glob
import logging
import os
from pathlib import Path
import shutil

import pandas as pd
import vcf  # pyvcf

from src import allele_frequency_filter, windowed_snp_counter
from src import zscore_distribution, vcf_cross_check
from src import plot, alter_original_vcf

logger = logging.getLogger(__name__)
def init_logger(log_filename):
    logger.setLevel(logging.INFO)
    file_formatter = logging.Formatter('%(asctime)s: %(message)s')
    stream_formatter = logging.Formatter('%(message)s')
    file_handler = logging.FileHandler(log_filename)
    stream_handler = logging.StreamHandler()
    file_handler.setFormatter(file_formatter)
    stream_handler.setFormatter(stream_formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger

# ---------------------- Helper Functions ---------------------
def output_directory_setup(output_directory, project, window_size_str):
    # Make log output directory
    Path(f"{output_directory}/{project}_{window_size_str}/logs").mkdir(parents=True, exist_ok=True)

    # Allele Frequency Filter output
    Path(f"{output_directory}/{project}_{window_size_str}/candidate_sites").mkdir(parents=True, exist_ok=True)

    # VCF_Check output
    Path(f"{output_directory}/{project}_{window_size_str}/dropped_snp_info").mkdir(parents=True, exist_ok=True)

    # Windowed SNP counter output
    Path(f"{output_directory}/{project}_{window_size_str}/windowed_snp_counts").mkdir(parents=True, exist_ok=True)

    # Z-score dfistribution outputs
    Path(f"{output_directory}/{project}_{window_size_str}/zscore_distribution").mkdir(parents=True, exist_ok=True)
    Path(f"{output_directory}/{project}_{window_size_str}/zscore_boolean").mkdir(parents=True, exist_ok=True)

    # Make project figures output directories
    Path(f"{output_directory}/{project}_{window_size_str}/figures/z_score_distributions/").mkdir(parents=True, exist_ok=True)
    Path(f"{output_directory}/{project}_{window_size_str}/figures/z_score_manhattan_plots").mkdir(parents=True, exist_ok=True)
    Path(f"{output_directory}/{project}_{window_size_str}/figures/z_score_heatmaps/html").mkdir(parents=True, exist_ok=True)
    Path(f"{output_directory}/{project}_{window_size_str}/figures/z_score_heatmaps/svg").mkdir(parents=True, exist_ok=True)

    # Output stats
    Path(f"{output_directory}/{project}_{window_size_str}/stats/").mkdir(parents=True, exist_ok=True)
    return


def convert_window_size(ws):
    """
    This function converts the shorthand input window size
    and returns an integer of the same value (i.e. "100kb" == int(100000))
    Args:
        ws: window size (bp/kb/mb)
    Returns: Integer of window size
    """
    window_size = None
    if "bp" in ws:
        window_size = int(ws.strip("bp"))*100
    elif "kb" in ws:
        window_size = int(ws.strip("kb"))*1000
    elif "mb" in ws:
        window_size = int(ws.strip("mb"))*10000
    return window_size


def listdir_nohidden(path):
    """This function takes in a pathway to a directory
       and collects all files expect hidden files."""
    return glob.glob(os.path.join(path, '*'))


def sample_ref_dict_maker(samp_ref_df):
    samp_ref_dict = {}
    for sample, ref_num, pop1 in zip(
            samp_ref_df['Sample'],
            samp_ref_df['Subspecies_Reference_Number'],
            samp_ref_df['PG1']
    ):
        samp_ref_dict[str(sample)] = (ref_num, pop1)
        continue
    return samp_ref_dict


def init_iid_log(log_filename):
    if log_filename.exists():
        return
    else:
        logger.info("=======================================")
        logger.info(" --- Initializing Introgression ID --- ")
        logger.info("=======================================")
        return

# ----------------------- Main Function ----------------------
def introgression_id():
    parser = argparse.ArgumentParser(
        description="This tool is used for the bison-cattle introgression study.",
    )
    # --- Stage Arguments ---
    parser.add_argument(
        '--allele',
        action="store_true",
        help='Only run allele frequency step',
        default=False,
    )
    parser.add_argument(
        '--vcf_check',
        action="store_true",
        help='Runs introgressed positions against vcf file of cattle variant data.',
        default=False,
    )
    parser.add_argument(
        '--count',
        action="store_true",
        help='Only run SNP counting step',
        default=False,
    )
    parser.add_argument(
        '--zscore',
        action="store_true",
        help='Only run standardization step',
        default=False,
    )
    parser.add_argument(
        '--plot',
        action="store_true",
        help='Plot data',
        default=False,
    )
    parser.add_argument(
        '--alter',
        action="store_true",
        help='Update original vcf files with introgression information',
        default=False,
    )
    # --- Required Data Arguments ---
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        action='store',
        required=True,
        help='Indicate which chromosome to run.',
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        action='store',
        help='Output directory pathway',
        default="output_data"
    )
    parser.add_argument(
        '-p',
        '--project',
        type=str,
        action='store',
        required=True,
        help='Indicate a project name. This will keep runs with different data separate.',
    )
    parser.add_argument(
        '-r',
        '--reference',
        type=str,
        action='store',
        default=None,
        help='Indicate pathway to standard reference excel sheet.',
    )
    parser.add_argument(
        '-b',
        '--bed',
        type=str,
        action='store',
        default=None,
        help='Indicate which filter to use.',
    )
    parser.add_argument(
        '-w',
        '--window_size',
        type=str,
        action='store',
        help='Indicate which filter to use.',
    )
    parser.add_argument(
        '-t',
        '--threshold',
        type=float,
        action='store',
        help='Z-score threshold.',
    )
    parser.add_argument(
        '-v',
        '--vcf_file',
        type=str,
        action='store',
        help='Cattle VCF file',
    )
    
    # -- Optional Arguments --
    parser.add_argument(
        '--auto_graph',
        action="store_true",
        help="Automatically open all graphs in default browser (NOTE: If you have a large number of samples, I would advise against using this feature as it will open a graph for every sample provided)",
        default=False,
    )
    parser.add_argument(
        '--vcf_chrom',
        type=str,
        action='store',
        default=None,
        help='Chromosome number when VCF to crosscheck is not entire genome',
    )
    args = parser.parse_args()

    input_directory = args.input
    output_directory = args.output
    allele = args.allele
    vcf_check = args.vcf_check
    count = args.count
    zscore = args.zscore
    plot_data = args.plot
    alter_vcf = args.alter
    project = args.project
    sample_reference_info = args.reference
    chromosome_length_bed_file = args.bed
    window_size_str = args.window_size  # 10bp/kb/mb
    zscore_threshold = args.threshold
    auto_graph = args.auto_graph
    cattle_vcf_file = args.vcf_file
    vcf_chrom = args.vcf_chrom

    # --- Convert window size to int ---
    window_size_int = convert_window_size(window_size_str)

    # --- Set up directories ---
    output_directory_setup(output_directory, project, window_size_str)
    log_output_dir = Path(f"{output_directory}/{project}_{window_size_str}/logs")
    log_filename = log_output_dir / f"{project}_{window_size_str}.log"
    init_iid_log(log_filename)
    logger = init_logger(log_filename)

    # Dict organization={"sample": {"LIR": [], "SIR": [], "NIS": VALUE}}
    sample_stats_keeper = {}

    # If no argument specified, run all functions
    args_list = [allele, count, zscore, vcf_check, plot_data, alter_vcf]

    if True not in args_list:
        allele = True
        count = True
        zscore = True
        vcf_check = True
        plot_data = True
        alter_vcf=True
        print()
        print("** No process specified, running entire pipeline")
        print()
        pass
    # Filter data to retain potentially introgressed positions
    if allele:
        print()
        print(f"| --- Running allele frequency filter --- ")
        print()
        # --- Stats tracking ---
        vcf_file_stats = {
            "Chromosome": [],
            "Total_Records": [],
            "Number_of_Introgressed_SNPs": [],
        }
        # Colelct chromosome VCF files
        chromosome_files = [Path(f) for f in listdir_nohidden(input_directory)]
        for chrom_file in chromosome_files:
            print()
            print(f"|---- Running AAF for {chrom_file.name} ---- ")
            # Log which chromosome is being run currently
            logger.info(f"--- Running {chrom_file.name} ---")
            # Pull chromosome name from file
            chromosome_number = chrom_file.stem.lower()
            # load in vcf file
            vcf_reader = vcf.Reader(filename=str(chrom_file))
            # load in sample reference data
            samp_ref_df = pd.read_excel(sample_reference_info)
            number_of_bison = len(samp_ref_df[samp_ref_df["Subspecies_Reference_Number"] == 1])
            number_of_cattle = len(samp_ref_df[samp_ref_df["Subspecies_Reference_Number"] == 0])
            # Create sample reference dictionary for speedy look-up.
            samp_ref_dict = sample_ref_dict_maker(samp_ref_df)
            # If True, run Allele Frequency Filter

            init_count_of_records = 0
            introgressed_record_total = 0

            output_location = f"{output_directory}/{project}_{window_size_str}/candidate_sites/" \
                            f"{chromosome_number}.introgressed_snps.tsv"

            allele_freq_dict = {
                'Chromosome': [],
                'Position': [],
                'Reference_Base': [],
                'Alternative_Base': [],
                'Reference_Frequency': [],
                'Alternative_Frequency': [],
                'Call_Rate': [],
                'cow_hom_ref': [],
                'bison_hom_ref': [],
                'wb_hom_ref': [],
                'cow_hom_alt': [],
                'bison_hom_alt': [],
                'wb_hom_alt': [],
                'cow_het': [],
                'bison_het': [],
                'wb_het': [],
                'total_cows_called': [],
                'total_bison_called': [],
                'het_bison_samples': [],
                'PG1': [],
            }
            # Each record is checked against filter criteria
            # If it passes it's output is added to the a collection dictionary
            for record in vcf_reader:
                init_count_of_records += 1
                result = allele_frequency_filter.allele_frequency_filter(
                    record,
                    samp_ref_dict,
                    number_of_bison,
                    number_of_cattle,
                )
                if result is None:
                    continue
                else:
                    introgressed_record_total += 1
                    for key, value in zip(allele_freq_dict.keys(), result):
                        allele_freq_dict[key].append(value)
                        continue
            # Load allele_frequency_dict into a DataFrame and output to tsv file
            allele_freq_df = pd.DataFrame.from_dict(allele_freq_dict)
            allele_freq_df.to_csv(output_location, index=False, sep='\t', mode='w')
            # Return dictionary of collected stats
            output_dict = {
                "init_count_of_records": init_count_of_records,
                "introgressed_record_total": introgressed_record_total,
            }

            # -------- Record Statistics -------
            vcf_file_stats["Chromosome"].append(chromosome_number)
            vcf_file_stats["Total_Records"].append(output_dict["init_count_of_records"])
            vcf_file_stats["Number_of_Introgressed_SNPs"].append(output_dict["introgressed_record_total"])
            logger.info(f"Total number of records in VCF file {chrom_file.name}= {output_dict}")
            print(f"|---- Total number of records in VCF file {chrom_file.name}= {output_dict}")
            print()
        vcf_output_location = f"{output_directory}/{project}_{window_size_str}/candidate_sites/" \
                            f"candidate_site_counts.tsv"
        vcf_stats_df = pd.DataFrame(data=vcf_output_location)
        vcf_stats_df.to_csv(vcf_output_location, sep='\t')
        pass
    # Cross check filter1 output to additional variant file
    if vcf_check:
        cattle_vcf_file = Path(cattle_vcf_file)
        logger.info(f"| ---- Secondary VCF Validation: {cattle_vcf_file.name} ---- |")
        dropped_snp_output_dir = Path(f"{output_directory}/{project}_{window_size_str}/dropped_snp_info/")
        candidate_positions_path = Path(f"{output_directory}/{project}_{window_size_str}/candidate_sites/")
        candidate_positions_files = sorted([f for f in Path(candidate_positions_path).iterdir() if f.is_file()])
        
        if vcf_chrom:
            pass
        elif len(candidate_positions_files) > 1:
            old_positions_files = [f for f in Path(candidate_positions_path).iterdir() if f.is_file()]
            now = datetime.now()
            cdate = now.strftime("%m_%d_%y")
            old_data_output_dir = candidate_positions_path / f"previous_candidate_sites_{cdate}/"
            old_data_output_dir.mkdir(parents=True, exist_ok=True)
            for f in old_positions_files:
                new_file_name = old_data_output_dir / f.name
                shutil.move(f.as_posix(), new_file_name.as_posix())
            candidate_positions_files = sorted([f for f in old_data_output_dir.iterdir() if f.is_file()])
        introgressed_output_dir = Path(f"{output_directory}/{project}_{window_size_str}/candidate_sites/")        
        # load in vcf file
        vcf_reader = vcf.Reader(filename=cattle_vcf_file.as_posix())
        vcf_cross_check.vcf_check(
            candidate_positions_files,
            vcf_reader,
            introgressed_output_dir,
            dropped_snp_output_dir,
            vcf_chrom,
            logger,
        )
        pass
    # Count number of SNPs per window
    if count:
        print()
        print("|---- Running Windowed SNP Counter ----")
        print()
        # Per-sample dictionary of DataFrames [key=sample: value=pd.DataFrame()]
        sample_snp_counts = dict()
        bison_sample_stats = dict()
        # Collect filtered files into a list
        file_path=Path(f"{output_directory}/{project}_{window_size_str}/candidate_sites/")
        files = [f for f in file_path.iterdir() if f.is_file()]
        # Iterate through each file (chromosome) and feed it to calc_window_denisty
        for current_file in files:
            sliding_window_returns = windowed_snp_counter.windowed_snp_counter(
                sample_snp_counts,
                current_file,
                chromosome_length_bed_file,
                window_size_int,
            )
            sample_snp_counts = sliding_window_returns[0]
            # Update sample stats
            for sample in sliding_window_returns[1].keys():
                try:
                    bison_sample_stats[sample] += sliding_window_returns[1][sample]
                except KeyError:
                    bison_sample_stats[sample] = sliding_window_returns[1][sample]

        # Output completed dataframe to tsv file
        for bison_sample in sample_snp_counts.keys():
            output_file = f"{output_directory}/{project}_{window_size_str}/windowed_snp_counts/" \
                        f"{bison_sample}.filter1.filter2.snp_counts.window_size_{window_size_str}.tsv"
            sample_snp_counts[bison_sample].to_csv(output_file, sep="\t")

        # Record sample stats
        for bison_sample in bison_sample_stats.keys():
            sample_stats_keeper[bison_sample] = {
                "Longest Introgressed Region": None,
                "Shortest INtrogressed Region": None,
                "Number of Introgressed Positions": bison_sample_stats[bison_sample],
            }
        pass
    # Z-score tranform data
    if zscore:
        print()
        print("|---- Running Z-score Distribution ---- ")
        print()
        snp_count_dist_stats = zscore_distribution.zscore_distribution(
            output_directory,
            project,
            chromosome_length_bed_file,
            window_size_int,
            window_size_str,
            zscore_threshold,
            auto_graph,
        )

        stats_df = pd.DataFrame(
            columns=[
                "Sample",
                "Number of introgressed bases",
                "Number of introgressed windows",
                "Longest introgressed region",
                "Shortest introgressed region",
                "Percent introgressed per base",
                "Percent introgressed per window",
            ]
        )

        for index, sample in enumerate(snp_count_dist_stats.keys()):
            stats_df.at[index, "Sample"] = sample
            stats_df.at[index, "Number of introgressed bases"] = snp_count_dist_stats[sample]["number_of_introgressed_bases"]
            stats_df.at[index, "Number of introgressed windows"] = snp_count_dist_stats[sample]["number_of_introgressed_windows"]
            stats_df.at[index, "Longest introgressed region"] = snp_count_dist_stats[sample]["LIR"]
            stats_df.at[index, "Shortest introgressed region"] = snp_count_dist_stats[sample]["SIR"]
            stats_df.at[index, "Percent introgressed per base"] = snp_count_dist_stats[sample]["percent_introgressed_per_base"]
            stats_df.at[index, "Percent introgressed per window"] = snp_count_dist_stats[sample]["percent_introgressed_per_window"]
            continue

        stats_df.to_csv(
            f"{output_directory}/{project}_{window_size_str}/stats/per_sample_stats_{window_size_str}.tsv",
            sep="\t",
            index=False,
        )
        pass
    # Create per-sample heatmaps
    if plot_data:
        print()
        print("|---- Creating per-sample introgression heatmaps ---- ")
        plot.plot_data(
            output_directory,
            project,
            window_size_str,
            window_size_int,
            zscore_threshold,
            auto_graph,
        )

    # Add introgression info to original vcf files
    if alter_vcf:
        # File pathways
        init_vcf_path = Path(input_directory)
        # Collect files
        original_vcfs = sorted([f for f in init_vcf_path.iterdir() if f.is_file()])
        # Run each original chromosome file
        for chrom_file in original_vcfs:
            chromosome_name = chrom_file.stem
            logger.info("")
            logger.info(f"|---- Updating original VCF files with introgression data for chromosome {chromosome_name} ---- ")
            logger.info("")
            alter_original_vcf.alter_original_vcf(
                input_directory,
                output_directory,
                project,
                window_size_str,
                window_size_int,
                chrom_file,
                chromosome_name,
                sample_reference_info,
                logger,
            )
        return

    return


if __name__ == '__main__':
    introgression_id()

