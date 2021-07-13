import logging
from pathlib import Path

import numpy as np
import pandas as pd
import plotly
import plotly.express as px
from scipy.stats import zscore, norm

# ---------------------- Plotting Functions ----------------------

def plot_normpdf(
    z_score_list,
    output_directory,
    project,
    window_size_str,
    auto_graph=False,
    title_indicator=None,
    all_samples=False,
):
    """Plot normal probability density curve of z-score transformed data"""
    z_score_list.sort()
    z_mean = np.mean(z_score_list)
    z_std = np.std(z_score_list)
    pdf = norm.pdf(z_score_list, loc=z_mean, scale=z_std)
    fig = px.line(
        x=z_score_list,
        y=pdf,
    )
    # Make output directory
    if all_samples:
        figure_filename = f"{output_directory}/{project}_{window_size_str}/figures/z_score_distributions/all_samples_pdf.html"
    else:
        figure_filename = f"{output_directory}/{project}_{window_size_str}/figures/z_score_distributions/{title_indicator}.html"
    plotly.offline.plot(
        fig,
        filename=figure_filename,
        auto_open=auto_graph
    )
    return

# ----------------------- Helper Functions ------------------------
def set_nan(df, chrom_bed_file):
    """This function will take in a dataframe and chromosome length bed file
    and will replace 0's with np.nan according to each chromosome length.
    This will fix any issues when calculating Z-scores"""
    # Build dictionary of key=chromosome and value=chromosome_length
    chrom_length_dict = {}
    for v in chrom_bed_file.itertuples():
        chrom_length_dict[v[1]] = v[2]
        continue
    # Iterate through each column
    for chrom in df.columns.to_list():
        current_chrom_length = chrom_length_dict[str(chrom)]

        # Iterate through each value of a column in reverse
        for index, value in zip(
                reversed(df.index.to_list()),
                reversed(df[chrom].to_list())
        ):
            # Check if index is greater than length of chromosome
            if index > current_chrom_length:
                df.at[index, chrom] = np.nan
            else:
                break
    return df


def dataframe_to_list(dataframe):
    """
    This function takes in a dataframe of zscore values and puts
    all values into a single list for plotting. It omits all NaN values
    since they represent non-existing positions on chromosomes smaller than
    the largest chromosome.
    @param dataframe: Zscore dataframe to be put into list
    @return: List of zscore values for entire genome
    """
    zscore_list = []
    for col in dataframe.columns.to_list():
        for value in dataframe[col]:
            if str(value) == str(np.nan):  # If value is NaN, skip value
                continue
            else:
                zscore_list.append(value)
    return zscore_list


def calc_mean(list_of_values):
    """Calculate mean of list"""
    return sum(list_of_values) / len(list_of_values)


def calc_stdev(list_of_values):
    """Calculate standard deviation"""
    return np.std(list_of_values)


def make_zscore_df(dataframe, pop_mean, stdev):
    """This Function z-score transforms the provided dataframe"""
    zscore_df = dataframe.applymap(lambda x: (x - pop_mean)/stdev)
    return zscore_df


def update_all_sample_df(sample_snp_counts_df, all_sample_df):
    """This function iterates through sample snp count df
    and updates all_sample_df with it's snp count data"""
    for col in sample_snp_counts_df.columns.to_list():
        for index in sample_snp_counts_df.index.to_list():
            sample_window_value = sample_snp_counts_df.at[index, col]
            if str(sample_window_value) == str(np.nan):
                continue
            else:
                all_sample_df.at[index, col] += float(sample_window_value)
                continue
    return all_sample_df


def check_window_is_above_threshold(value, zscore_threshold):
    """Check if provided value is above z-score threshold"""
    if value > zscore_threshold:
        return True
    else:
        return False


def count_introgressed_positions(df, zscore_boolean_df):
    """Iterate through sample windowed SNP count df and introgression boolean
    df and tally the number of SNP's in windows with a zscore above the given
    threshold."""
    total_per_base = 0
    total_per_window = 0
    for col in df.columns.to_list():
        for v1, v2 in zip(df[col], zscore_boolean_df[col]):
            if v2:
                total_per_window += 1
                total_per_base += int(v1)
                continue
            else:
                continue
    return total_per_base, total_per_window


def find_longest_introgressed_region(zscore_boolean_df, window_size):
    """Identifies the longest continuous run of True windows and returns
    the the chromomose, start, stop, and length of the introgressed region"""

    longest_chrom = None
    longest_start = 0
    longest_stop = 0

    pos_flag = False
    pos_chrom = None
    pos_start = 0
    pos_stop = 0

    df_cols = zscore_boolean_df.columns.to_list()
    df_index = zscore_boolean_df.index.to_list()

    # Iterate through each chromosome
    for chrom in df_cols:
        # Run through the windows (index) and boolean value
        for index, value in zip(df_index, zscore_boolean_df[chrom]):
            if value:
                if not pos_flag:
                    pos_flag = True
                    pos_chrom = str(chrom)
                    pos_start = int(index)
                    pos_stop = int(index)
                    continue
                elif pos_flag:
                    pos_stop = int(index)
            else:
                pos_diff = (pos_stop - pos_start)
                longest_diff = (longest_stop - longest_start)
                if pos_flag:
                    if pos_diff < longest_diff:
                        pos_chrom = None
                        pos_start = 0
                        pos_stop = 0
                        pos_flag = False
                        continue
                    elif pos_diff > longest_diff:
                        longest_chrom = pos_chrom
                        longest_start = pos_start
                        longest_stop = pos_stop
                        pos_chrom = None
                        pos_start = 0
                        pos_stop = 0
                        pos_flag = False
                        continue
                    elif (pos_diff == 0) and (pos_start > 0) and (pos_stop > 0):
                        longest_chrom = pos_chrom
                        longest_start = pos_start
                        longest_stop = pos_stop
                        pos_chrom = None
                        pos_start = 0
                        pos_stop = 0
                        pos_flag = False
                else:
                    # No current possiblle region
                    continue
    out = [
        (
            longest_chrom,
            longest_start,
            longest_stop,
            (longest_stop-longest_start+window_size)
        )
    ]
    return out


def find_shortest_region(zscore_boolean_df, window_size):
    """This function will run through the introgressed boolean df for
    a given bison sample. It finds the longest introgressed region
    and returns the chromosome, start position, stop position, and
    the length of the region."""
    shortest_chrom = None
    shortest_start = 0
    shortest_stop = 0

    pos_flag = False
    pos_chrom = None
    pos_start = 0
    pos_stop = 0

    df_cols = zscore_boolean_df.columns.to_list()
    df_index = zscore_boolean_df.index.to_list()

    # Iterate through each chromosome
    for chrom in df_cols:
        # Run through the windows (index) and boolean value
        for index, value in zip(df_index, zscore_boolean_df[chrom]):
            if value:
                if not pos_flag:
                    pos_flag = True
                    pos_chrom = str(chrom)
                    pos_start = int(index)
                    pos_stop = int(index)
                    continue
                elif pos_flag:
                    pos_stop = int(index)
            else:
                pos_diff = (pos_stop - pos_start)
                shortest_diff = (shortest_stop - shortest_start)
                if pos_flag:
                    if pos_diff > shortest_diff:
                        pos_chrom = None
                        pos_start = 0
                        pos_stop = 0
                        pos_flag = False
                        continue
                    elif pos_diff < shortest_diff:
                        shortest_chrom = pos_chrom
                        shortest_start = pos_start
                        shortest_stop = pos_stop
                        pos_chrom = None
                        pos_start = 0
                        pos_stop = 0
                        pos_flag = False
                        continue
                    elif (pos_diff == 0) and (pos_start > 0) and (pos_stop > 0):
                        shortest_chrom = pos_chrom
                        shortest_start = pos_start
                        shortest_stop = pos_stop
                        pos_chrom = None
                        pos_start = 0
                        pos_stop = 0
                        pos_flag = False
                else:
                    # No current possiblle region
                    continue
    out = [
        (
            shortest_chrom,
            shortest_start,
            shortest_stop,
            (shortest_stop-shortest_start+window_size)
        )
    ]
    return out

# ----------------------- Main Function Call ------------------------
def zscore_distribution(
        output_directory,
        project,
        chromosome_length_bed_file,
        window_size_int,
        window_size_str,
        zscore_threshold,
        auto_graph,
):
    """This function will z-score transform the windowed snp count files
    for each sample. It will also make a collective dataframe of all samples
    results to obtain a z-score cut-off to use as the significance threshold.
    """
    # Initiate zscore_distribution log file
    logging.basicConfig(filename=F'logs/{project}_{window_size_str}_zscore_distribution.log', level=logging.INFO, filemode='w', format='')
    logging.info('')
    logging.info("|---- Running Z-score Distribution ---- ")
    logging.info('')
    # Load chromosome bed file information into DataFrame
    read_chromosome_bed_file = pd.read_csv(
        chromosome_length_bed_file,
        sep="\t",
        thousands=",",
        comment="#",
        names=["chromosome", "length"],
        dtype={"chromosome": str, "length": int},
    )
    length_of_genome = sum([v for v in read_chromosome_bed_file["length"]])
    # Input pathways + files
    count_pathway = f"{output_directory}/{project}_{window_size_str}/windowed_snp_counts/"
    raw_files = [f for f in Path(count_pathway).iterdir() if (f.is_file()) and (f.stem[0] != ".")]
    # Output pathways
    zscore_output_filename = f"{output_directory}/{project}_{window_size_str}/zscore_distribution/" \
                             f"all_samples_cumulative_zscores.tsv"

    # Sample data collection DataFrame
    all_sample_df = pd.DataFrame()
    sample_dfs = dict()

    # Iterate through each sample file to collect population data
    for sample_file in raw_files:
        # Load files into dataframes for processing
        sample_snp_counts_df = pd.read_csv(
            sample_file,
            sep="\t",
            index_col=[0],
        )
        # Collect file metadata
        sample_name = str(sample_file.name.split(".")[0])
        # Add sample_df to sample_dfs
        sample_dfs[sample_name] = sample_snp_counts_df
        # Initiate all sample dataframe is empty
        if all_sample_df.empty:
            # Set first input data
            all_sample_df = sample_snp_counts_df
            continue
        else:
            # Update all_sample_df
            update_all_sample_df(sample_snp_counts_df, all_sample_df)
            continue
        continue

    # zscore transform population data and output to file
    df_list = dataframe_to_list(all_sample_df)
    population_mean = calc_mean(df_list)
    population_stdev = calc_stdev(df_list)

    print(f"Population Mean = {population_mean}")
    print(f"Population Stdev = {population_stdev}")
    print()
    logging.info(f"Population Mean = {population_mean}")
    logging.info(f"Population Stdev = {population_stdev}")
    logging.info('')

    all_sample_zscore_df = make_zscore_df(all_sample_df, population_mean, population_stdev)
    all_sample_zscore_df.to_csv(zscore_output_filename, sep="\t")
    # Plot normal probability density function
    plot_normpdf(
        dataframe_to_list(all_sample_zscore_df),  # Converts DataFrame into list w/o np.nan's
        output_directory,
        project,
        window_size_str,
        auto_graph=auto_graph,
        title_indicator="All Samples",
        all_samples=True,
    )
    # Collect Statistics
    sample_stats = {}
    # Go back through and make per-sample z-score df's
    for sample_name in sample_dfs.keys():
        sample_snp_counts_df = sample_dfs[sample_name]
        sample_df_list = dataframe_to_list(sample_snp_counts_df)
        sample_stdev = calc_stdev(sample_df_list)
        # Output file names
        zscore_output_filename = f"{output_directory}/{project}_{window_size_str}/zscore_distribution/{sample_name}.zscore.tsv"
        zb_output = f"{output_directory}/{project}_{window_size_str}/zscore_boolean/{sample_name}.zscore_boolean.tsv"
        # Set up sample key-value pair in stats collection
        if sample_stats.get(sample_name) == None:
            sample_stats[sample_name] = {
                "number_of_introgressed_bases": 0,
                "number_of_introgressed_windows": 0,
                "LIR": [],
                "SIR": [],
                "percent_introgressed_per_base": 0,
                "percent_introgressed_per_window": 0,
            }

        # Make zscore DataFrame + output
        zscore_df = make_zscore_df(sample_snp_counts_df, population_mean, sample_stdev)
        zscore_df.to_csv(zscore_output_filename, sep="\t")

        # Convert zscore df to True if above zscore threshold + output boolean df
        zscore_boolean_df = zscore_df.applymap(lambda x: check_window_is_above_threshold(x, zscore_threshold))
        zscore_boolean_df.to_csv(zb_output, sep="\t")

        # Update per-sample % introgressed
        number_of_introgressed_sites = count_introgressed_positions(sample_snp_counts_df, zscore_boolean_df)

        # Find the longest and shortest block of introgressed postions
        longest_introgressed_region = find_longest_introgressed_region(zscore_boolean_df, window_size_int)
        shortest_introgressed_region = find_shortest_region(zscore_boolean_df, window_size_int)

        # Add stats to output dictionary
        sample_stats[sample_name]["percent_introgressed_per_base"] = round(
            (number_of_introgressed_sites[0] / length_of_genome) * 100,
            8
        )
        sample_stats[sample_name]["percent_introgressed_per_window"] = round(
            (number_of_introgressed_sites[1] / (length_of_genome / window_size_int)) * 100,
            4
        )
        sample_stats[sample_name]["LIR"] += longest_introgressed_region
        sample_stats[sample_name]["SIR"] += shortest_introgressed_region
        sample_stats[sample_name]["number_of_introgressed_bases"] = int(number_of_introgressed_sites[0])
        sample_stats[sample_name]["number_of_introgressed_windows"] = int(number_of_introgressed_sites[1])

        print(f"Z-score transformed {sample_name}")

        continue
    print(f"Completed population data")
    logging.info(f"Completed population data")
    return sample_stats
