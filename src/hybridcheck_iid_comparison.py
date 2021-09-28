"""
Author: Andrew Harris
Python 3.8.3
This script will map the Hybrid check data where the percentage of BosTau_BisBis is larger than
the percentage of BisBis_BisBis on a heatmap
"""
import argparse
from pathlib import Path
import os

import numpy as np
import pandas as pd
import plotly
import plotly.graph_objects as go


def _make_outdirs(OUTPUT):
    hc_heatmaps = OUTPUT / "HC_output_heatmaps"
    hc_heatmaps.mkdir(parents=True, exist_ok=True)
    hc_iid_heatmaps = OUTPUT / "HC_IID_comparison_output_heatmaps"
    hc_iid_heatmaps.mkdir(parents=True, exist_ok=True)
    return hc_heatmaps, hc_iid_heatmaps


def round_up(x, window_size):
    """
    @status: PASSING

    This function takes in a position and rounds it up to the nearest 100,000th position.
    @param window_size: Window size used in analysis
    @param x: Position of record
    @return:Position rounded up to the nearest 100,000th.
    """
    return x if x % window_size == 0 else x + window_size - x % window_size


def set_nan(df, chrom_bed_file):
    """This function will take in a dataframe and chromosome length bed file
    and will replace 0's with np.nan according to each chromosome length.
    This will fix any issues when calculating Z-scores.

    Args:
        df ([Object]): [sample snp count dataframe]
        chrom_bed_file ([object]): [Pandas dataframe with chromosome length information]
    """
    # Build dictionary of key=chromosome and value=chromosome_length
    chrom_length_dict = {}
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
                df.at[index, chromosome] = np.nan
    return df


def add_plotting_column(df):
    """
    This function will create a new column ("plot_value") in the dataframe
    and will set either a 0 if BisBis_BisBis % > BosTau_BisBis %, and 1 if
    BisBis_BisBis % < BosTau_BisBis %
    Args:
        df: input file DataFrame

    Returns: new dataframe with added "plot_value" column
    """
    df["plot_value"] = [0]*len(df)
    for index, row in enumerate(df.itertuples()):
        BosTau_BisBis = row[3]
        BisBis_BisBis = row[4]

        if BosTau_BisBis > BisBis_BisBis:
            df["plot_value"] = 1.0
        else:
            continue
    return df


def extend_position_data(plot_df, window_size, longest_chromosome, num_of_chroms):
    cols = [c for c in range(1, (num_of_chroms + 1))]
    index = [n for n in range(window_size, (longest_chromosome + window_size), window_size)]
    new_df = pd.DataFrame(0.0, columns=cols, index=index)

    # Iterate through plotting df and update all values in roundup(start) roundup(stop) range in new_df
    for row in plot_df.itertuples():
        curr_chrom = int(row[1].strip("chr"))
        start_ru = round_up(int(row[2]), window_size)
        stop_ru = round_up(int(row[3]), window_size)
        plotting_value = row[6]

        for window in range(start_ru, (stop_ru + window_size), window_size):
            if plotting_value == 1.0:
                new_df.at[window, curr_chrom] = plotting_value
            else:
                continue
    return new_df


def make_plot_df(df):
    collection_dict = {"Chromosome": [], "Position": [], "Introgressed": []}
    for col in df.columns.to_list():
        for index in df.index.to_list():
            collection_dict["Chromosome"].append(col)
            collection_dict["Position"].append(index)
            collection_dict["Introgressed"].append(df.at[index, col])
            continue
    output_df = pd.DataFrame(data=collection_dict)
    return output_df


def plot_hc_heatmap(plot_df, sample_name, hc_heatmaps, hc_iid_heatmaps, auto_graph, just_hc=True):

    fig = go.Figure(go.Heatmap(
        x=plot_df["Position"],
        y=plot_df["Chromosome"],
        z=plot_df["Introgressed"],
        ygap=2.5,
        colorscale=['rgb(235, 234, 234)', 'rgb(255, 19, 0)'],
        zmin=0,
        zmax=1,
        colorbar=dict(
            tick0=0,
            dtick=1
        )
    ))
    fig.update_layout(
        title=f'{sample_name} HybridCheck Heatmap',
        yaxis_nticks=29,
        xaxis=dict(
            tickmode='linear',
            tick0=0,
            dtick=10000000,
        ),
        yaxis=dict(
            tickmode='linear',
            tick0=0,
            dtick=1
        ),
        height=900,
        width=1250,
        template="simple_white",
        font=dict(
            family="Arial",
            size=15,
        ),
    )
    fig.update_traces(showscale=False)
    fig.update_yaxes(
        title_text='Chromosome',
        ticks="outside",
        tickson="boundaries",
        ticklen=5,
        linecolor='black',
        mirror=True,
    )

    fig.update_xaxes(
        title_text='Position',
        showgrid=False,
        linecolor='black',
        mirror=True,
    )
    if just_hc:
        html_dir = hc_heatmaps / "html"
        html_dir.mkdir(parents=True, exist_ok=True)
        svg_dir = hc_heatmaps / "svg"
        svg_dir.mkdir(parents=True, exist_ok=True)
        html_filename = html_dir / f"{sample_name}_HybridCheck_heatmap.html"
        svg_filename = svg_dir / f"{sample_name}_HybridCheck_heatmap.svg"
    else:
        html_dir = hc_iid_heatmaps / "html"
        html_dir.mkdir(parents=True, exist_ok=True)
        svg_dir = hc_iid_heatmaps / "svg"
        svg_dir.mkdir(parents=True, exist_ok=True)
        html_filename = html_dir / f"{sample_name}_HybridCheck_regions_in_common_heatmap.html"
        svg_filename = svg_dir / f"{sample_name}_HybridCheck_regions_in_common_heatmap.svg"

    fig.write_image(svg_filename, format='svg', engine='kaleido')
    plotly.offline.plot(
        fig,
        filename=html_filename,
        auto_open=auto_graph,
    )
    return


def compare_outputs(introgressionid_df, hc_df, auto_graph, sample_name, hc_heatmaps, hc_iid_heatmaps):
    """
    This function will compare the two dataframes and identify + plot the regions where both approaches
    identified a region as introgressed.
    Args:
        introgressionid_df: IntrogressionID sample output Boolean df
        hc_df: HybridCheck sample plotting df
    Returns: None (outputs heatmap)
    """

    def _count_similar_regions(hc_df):
        """Count number of regions with 1.0"""
        total = 0
        for _, c in hc_df.iterrows():
            if c.Introgressed == 1.0:
                total += 1
                continue
            else:
                continue
        return total

    # Change HybridCheck results to only reflect regions 
    # found in both HybridCheck and IntrogressionID
    for row in hc_df.itertuples():
        index = row[0]
        chromosome = row[1]
        position = row[2]
        introgressed = row[3]

        # If window is introgressed in IntrogressionID output, continue
        if introgressed == 0:
            continue
        elif (introgressed == 1) and (introgressionid_df.at[position, str(chromosome)]):
            continue
        # If not, change the Introgressed value from 1 to 0.
        else:
            if str(hc_df.at[index, "Introgressed"]) == str(np.nan):
                continue
            else:
                hc_df.at[index, "Introgressed"] = 0
    
    # Calculate number of regions in common between 
    # HybridCheck and IntrogressionID results
    total_common_regions = _count_similar_regions(hc_df)
    print(f"{sample_name}: {total_common_regions}")
    
    # Plot results if not empty
    if hc_df.empty:
        return
    else:
        plot_hc_heatmap(hc_df, f"{sample_name} IntrogressionID", hc_heatmaps, hc_iid_heatmaps, auto_graph, just_hc=False)
    return


def plot_hybridcheck_heatmap():
    
    parser = argparse.ArgumentParser(
        description="This tool is used to plot HybridCheck introgressex regions on a heatmap. It also has the option"
                    "to compare the HybridCheck results to the IntrogressionID output Boolean file.",
    )
    
    # -- Optional Arguments --
    parser.add_argument(
        '--output',
        action="store",
        type=str,
        help='Output location for graphs',
    )
    parser.add_argument(
        '--hc_input',
        action="store",
        type=str,
        help='HybridCheck Input file',
    )
    parser.add_argument(
        '--window_size',
        action="store",
        type=int,
        help='Window size',
    )
    parser.add_argument(
        '--chrom_len',
        action="store",
        type=int,
        help='Chromosome length bed file',
    )
    parser.add_argument(
        '--compare_to_introgressionid',
        action="store_true",
        help='Run comparison to IntrogressionID output',
        default=False,
    )
    parser.add_argument(
        '--introgressionid_file',
        action="store",
        type=str,
        help='IntrogressionID file to compare',
    )
    parser.add_argument(
        '--auto_graph',
        action="store_true",
        default=False,
        help='IntrogressionID file to compare',
    )
    args = parser.parse_args()
    # -- Input variables --
    input_file = args.hc_input
    window_size = args.window_size
    chromosome_length_bed_file = args.chrom_len
    compare_files = args.compare_to_introgressionid
    file_to_compare = args.introgressionid_file
    auto_graph = args.auto_graph
    try:
        OUTPUT = Path(args.output)
    except FileNotFoundError:
        print("Output location cannot be found... check and rerun")
        exit(1)
    # -- Make output directories --
    hc_heatmaps, hc_iid_heatmaps = _make_outdirs(OUTPUT)
    # -- Set info --
    sample_name = Path(input_file).stem.split("_")[0]
    file_df = pd.read_csv(input_file, sep="\t")
    # -- Load in chromosome bed file information --
    read_chromosome_bed_file = pd.read_csv(
        chromosome_length_bed_file,
        sep="\t",
        thousands=",",
        comment="#",
        names=["chromosome", "length"],
        dtype={"chromosome": str, "length": int},
    )
    longest_chromosome = read_chromosome_bed_file["length"].max()
    number_of_chromosomes = len(read_chromosome_bed_file["chromosome"])
    added_plot_column = add_plotting_column(file_df)
    extended_df = extend_position_data(added_plot_column, window_size, longest_chromosome, number_of_chromosomes)
    nan_fixed_df = set_nan(extended_df, read_chromosome_bed_file)
    plot_df = make_plot_df(nan_fixed_df)
    # -- Compare HC and IID results if requested --
    if compare_files:
        read_introgressionid_file = pd.read_csv(file_to_compare, sep="\t", index_col=[0])
        compare_outputs(read_introgressionid_file, plot_df.copy(), auto_graph, sample_name, hc_heatmaps, hc_iid_heatmaps)
        pass
    # -- Plot heatmap for only HybridCheck data and output --
    plot_hc_heatmap(plot_df, sample_name, hc_heatmaps, hc_iid_heatmaps, auto_graph)
    return


if __name__ == '__main__':
    plot_hybridcheck_heatmap()
