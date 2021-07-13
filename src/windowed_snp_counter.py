import pandas as pd


def round_up(x, window_size):
    """round position up to the nearest multiple of the window size"""
    return x if x % window_size == 0 else x + window_size - x % window_size


def set_nan(df, chrom_bed_file):
    """Set DataFrame windows above the length of a given chromoosome to nan"""
    chrom_length_dict = {}
    # Build dictionary of key=chromosome and value=chromosome_length
    for v in chrom_bed_file.itertuples():
        chrom_length_dict[v[1]] = v[2]
        continue

    # Iterate through each column
    for chromosome in df.columns.to_list():
        current_chrom_length = chrom_length_dict[str(chromosome)]
        # Iterate through each value of a column in reverse
        for index, value in zip(
            reversed(df.index.to_list()),
            reversed(df[chromosome].to_list())
        ):
            # Check if index is greater than length of chromosome
            if index > current_chrom_length:
                df.at[index, chromosome] = pd.NA
            else:
                continue
    return df


def windowed_snp_counter(sample_snp_counts, introgression_pos_file, chromosome_length_bed_file, window_size):
    """This function uses a sliding window approach to count the number of SNP per window"""
    # Load files into dataframes for processing
    read_hybrid_pos_file = pd.read_csv(
        introgression_pos_file,
        sep="\t",
        usecols=['Chromosome', 'Position', 'het_bison_samples'],
    )
    # Load in chromosome bed file information
    read_chromosome_bed_file = pd.read_csv(
        chromosome_length_bed_file,
        sep="\t",
        thousands=",",
        comment="#",
        names=["chromosome", "length"],
        dtype={"chromosome": str, "length": int},
    )

    # Collect bison sample names found in introgressed positions file
    bison_with_introgressed_pos = read_hybrid_pos_file["het_bison_samples"].unique()

    # Collect sample stats and return (dictionary: key=sample, value=Number of introgressed snps)
    stats_collection = {}

    # Create collection DataFrame in sample_snp_counts if one is not already made
    for sample in bison_with_introgressed_pos:
        try:
            sample_snp_counts[sample]  # Check if entry exist
        except KeyError:
            # Need to add additional window_size to largest window size becuase range stops one position short.
            max_window = round_up(read_chromosome_bed_file["length"].max(), window_size) + window_size
            bare_df = pd.DataFrame(
                data=0, 
                columns=[n for n in read_chromosome_bed_file["chromosome"]],
                index=[i for i in range(window_size, max_window, window_size)],
                dtype=object
            )
            nan_fixed_df = set_nan(bare_df, read_chromosome_bed_file)
            sample_snp_counts[sample] = nan_fixed_df

    # Iterate through each unique sample and tally up count in n-(bp/kb/mb) windows.
    for bison_sample in bison_with_introgressed_pos:
        bison_introgressed_positions = read_hybrid_pos_file[read_hybrid_pos_file["het_bison_samples"] == bison_sample]
        # Record number of introgressed snps for current sample
        stats_collection[bison_sample] = len(bison_introgressed_positions)
        # Iterate through introgressed positions and add to rolling tally in dataframe
        for chromosome, pos, sample in zip(
            bison_introgressed_positions["Chromosome"],
            bison_introgressed_positions["Position"],
            bison_introgressed_positions["het_bison_samples"],
        ):
            sample_snp_counts[sample].at[round_up(pos, window_size), str(chromosome)] += 1
            continue
        continue
    return sample_snp_counts, stats_collection
