import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly
import plotly.express as px
import plotly.graph_objects as go
import plotly.offline as po
import plotly.io as pio

logger = logging.getLogger(__name__)
def init_logger(log_filename):
    logger.setLevel(logging.INFO)
    file_formatter = logging.Formatter('%(message)s')
    stream_formatter = logging.Formatter('%(message)s')
    file_handler = logging.FileHandler(log_filename)
    stream_handler = logging.StreamHandler()
    file_handler.setFormatter(file_formatter)
    stream_handler.setFormatter(stream_formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


# ------- Plotting Functions -------
def make_per_sample_heatmap(
    zscore_df,
    window_size_str,
    window_size_int,
    zscore_threshold,
    sample_name,
    output_directory,
    project,
    auto_graph,
):
    def make_plot_df(df):
        """Melt dataframe"""
        df.reset_index(level=0, inplace=True)
        melted_df = df.melt(id_vars=["index"], var_name="Chromosome", value_name="z-score", ignore_index=False)
        melted_df.columns = ["Position", "Chromosome", "z-score"]
        return melted_df


    def get_max_record(df):
        for _, record in df.iterrows():
            pos = record[0]
            chrom = record[1]
            z_score = record[2]
            output = f"Chromosome: {chrom}\nPosition: {pos}\nZ-score: {z_score}\n"
            return output

    def convert_to_presence(x, zscore_threshold):
        if x >= zscore_threshold:
            return 1
        else:
            return 0


    

    plot_df = make_plot_df(zscore_df)
    plot_df.dropna(inplace=True)
    above_threshold = plot_df[plot_df["z-score"] > zscore_threshold]
    if len(above_threshold) == 0:
        logger.info(f"Passing {sample_name}, no windows above min z-score threshold.")

    max_df_record = plot_df[plot_df["z-score"] == plot_df["z-score"].max()]
    min_df_record = plot_df[plot_df["z-score"] == plot_df["z-score"].min()]
    max_record = get_max_record(max_df_record)
    min_record = get_max_record(min_df_record)
    zscore_presence = plot_df["z-score"].apply(lambda x: convert_to_presence(x, zscore_threshold))

    # --- Binary Graph - Viridis ---
    fig = go.Figure(go.Heatmap(
        x=plot_df["Position"],
        y=plot_df["Chromosome"],
        z=zscore_presence,
        ygap=2.5,
        # colorscale=[[0, 'rgb(235, 234, 234)'], [1, 'rgb(255, 19, 0)']],
        colorscale="viridis",
        colorbar=dict(dtick=1, title='Presence'),
        zmin=0,
        zmax=1,
    ))
    # --- Binary Graph - Red/Grey ---
    fig2 = go.Figure(go.Heatmap(
        x=plot_df["Position"],
        y=plot_df["Chromosome"],
        z=zscore_presence,
        ygap=2.5,
        # colorscale=[[0, 'rgb(235, 234, 234)'], [1, 'rgb(255, 19, 0)']],
        colorscale=['rgb(235, 234, 234)', 'rgb(255, 19, 0)'],
        colorbar=dict(dtick=1, title='Presence'),
        zmin=0,
        zmax=1,
    ))

    # ---- Relative frequency plot ----
    # fig = go.Figure(go.Heatmap(
    #     x=plot_df["Position"],
    #     y=plot_df["Chromosome"],
    #     z=plot_df["z-score"],
    #     ygap=2.5,
    #     colorscale=[[0, 'rgb(235, 234, 234)'], [0.25, 'rgb(243, 143, 135)'], [0.5, 'rgb(241, 103, 92)'], [0.75, 'rgb(243, 70, 56)'], [1, 'rgb(255, 19, 0)']],
    #     zmin=5,
    #     zmax=max(plot_df["z-score"]),
    #     colorbar=dict(title='Z-score')
    # ))
    fig.update_traces(showscale=False)
    
    fig.update_layout(
        title=f'{sample_name} IntrogressionID Heatmap',
        xaxis=dict(
            title="Position",
            tickmode='linear',
            tick0=0,
            dtick=10000000,
            range=[0, 160000000]
        ),
        yaxis=dict(
            title="Chromosome",
            tickmode='linear',
            tick0=0,
            dtick=1
        ),
        template='simple_white',
        font=dict(
            family="Arial",
            size=15,
        ),
        height=900,
        width=1250,
    )
    fig2.update_layout(
        title=f'{sample_name} IntrogressionID Heatmap',
        xaxis=dict(
            title="Position",
            tickmode='linear',
            tick0=0,
            dtick=10000000,
            range=[0, 160000000]
        ),
        yaxis=dict(
            title="Chromosome",
            tickmode='linear',
            tick0=0,
            dtick=1
        ),
        template='simple_white',
        font=dict(
            family="Arial",
            size=15,
        ),
        height=900,
        width=1250,
    )
    fig2.update_traces(showscale=False)
    # HTML and SVG output file names - Viridis
    svg_filename_viridis = f"{output_directory}/{project}_{window_size_str}/figures/z_score_heatmaps/" \
                      f"svg/{sample_name}_heatmap_viridis.svg"
    html_filename_viridis = f"{output_directory}/{project}_{window_size_str}/figures/z_score_heatmaps/" \
                      f"html/{sample_name}_heatmap_viridis.html"
    # Write HTML and SVG images to files
    pio.write_image(fig=fig, file=svg_filename_viridis, format='svg', engine="kaleido")
    po.plot(
        fig,
        filename=html_filename_viridis,
        auto_open=auto_graph,
    )
    # HTML and SVG output file names - Red/Grey
    svg_filename_red_grey = f"{output_directory}/{project}_{window_size_str}/figures/z_score_heatmaps/" \
                      f"svg/{sample_name}_heatmap_red-grey.svg"
    html_filename_red_grey = f"{output_directory}/{project}_{window_size_str}/figures/z_score_heatmaps/" \
                      f"html/{sample_name}_heatmap_red-grey.html"
    # Write HTML and SVG images to files
    pio.write_image(fig=fig2, file=svg_filename_red_grey, format='svg', engine="kaleido")
    po.plot(
        fig2,
        filename=html_filename_red_grey,
        auto_open=auto_graph,
    )
    logger.info("")
    logger.info("==========================================")
    logger.info(f"--- {sample_name} ---")
    logger.info("")
    logger.info(f"Max z-score record: \n{max_record}")
    logger.info(f"Min z-score record: \n{min_record}")
    logger.info(f"Windows above z-score threshold: {len(above_threshold)}")
    logger.info("")
    logger.info("==========================================")
    return


def make_per_sample_manhattan_plot(
    sample_df,
    window_size_str,
    window_size_int,
    zscore_threshold,
    sample_name,
    output_directory,
    project,
    auto_graph,
):
    # Set up df datatypes
    melted_df = sample_df.melt(id_vars=["index"], var_name="Chromosome", value_name="z_score", ignore_index=False)
    melted_df.columns = ["Position", "Chromosome", "z_score"]
    melted_df.Chromosome = melted_df.Chromosome.astype('category')

    # Collect unique chromosomes
    unique_chroms = [str(c) for c in melted_df['Chromosome'].unique()]

    # Set chromosome categories
    melted_df.Chromosome = melted_df.Chromosome.cat.set_categories(unique_chroms, ordered=True)
    melted_df = melted_df.sort_values('Chromosome')
    melted_df.dropna(inplace=True)

    # Group data by chromosomes
    melted_df['ind'] = range(len(melted_df))
    melted_df_grouped = melted_df.groupby(by=['Chromosome'])

    fig = plt.figure(dpi=200, figsize=(15, 5))
    ax = fig.add_subplot(111)
    colors = ['red','blue']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(melted_df_grouped):
        group.plot(kind='scatter', x='ind', y='z_score',color=colors[num % len(colors)], ax=ax, s=3)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))

    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(melted_df)])
    ax.hlines(y=zscore_threshold, xmin=0, xmax=len(melted_df))
    ax.set_xlabel('Chromosome')
    manhattan_filename = f"./{output_directory}/{project}_{window_size_str}/figures/z_score_manhattan_plots/{sample_name}_manhattan.pdf"
    fig.savefig(manhattan_filename, dpi=300, orientation='landscape', format='pdf')
    plt.close()
    return


# ------- Main Function -------
def plot_data(output_directory, project, window_size_str,
              window_size_int, zscore_threshold, auto_graph):
    # Initiate plot log file
    filename=Path(f'{output_directory}/{project}_{window_size_str}/logs/{project}_{window_size_str}_plot.log')
    init_logger(filename)
    # Collect sample z-score files
    zscore_output_pathway = Path(f"{output_directory}/{project}_{window_size_str}/zscore_distribution/")
    sample_zscore_files = [f for f in zscore_output_pathway.iterdir() if (f.is_file()) and (f.stem[0] != ".")]
    # Create plots for sample
    for sample_file in sample_zscore_files:
        sample_name = sample_file.stem.split(".")[0]
        sample_df = pd.read_csv(sample_file, sep="\t", index_col=[0])
        make_per_sample_heatmap(
            sample_df,
            window_size_str,
            window_size_int,
            zscore_threshold,
            sample_name,
            output_directory,
            project,
            auto_graph,
        )
        make_per_sample_manhattan_plot(
            sample_df,
            window_size_str,
            window_size_int,
            zscore_threshold,
            sample_name,
            output_directory,
            project,
            auto_graph,
        )
        continue
    return
