#!/usr/bin/env python3

'''
This script contains a FUNCTION TO PLOT -log10 transformed p-values AGAINST log2FOLD CHANGES (FROM TRIPLICATE MEASUREMENTS) OF TWO DATASETS AS VOLCANO PLOT. The datasets need to be in the same input file. Therefore, the input file must be generated from a combination of following steps (in order): 
    1. format_proteomefile.py
    2. annotate.py (Only if annotations were missing in the original file)
    3. statistics.py

The input file needs to have 16 columns, where the extension '_1' is condition 1, and '_2' is condition 2: 
Accession_Number	Annotation_1	1_1 2_1 3_1 Annotation_2 1_2 2_2 3_2	Row_Average_1	STD_1	Row_Average_2	STD_2	Log2_Fold_Change	Adj_P_Value Log10_adj_pvalue


title = Title of the plot
x_axis_title = "log2 fold change"
y_axis_title = "-log10 pvalue"
point_radius = 4, can be changed but 4 works well
log2_fold_change = log2 fold change calculated from statistics.py function and saved in Log2_Fold_Change dataframe
transformed_pvalues = p-value calculated from statistics.py function and saved in Transformed_P_Value dataframe
annotations = annotations for the respective log2fold change and p-values. 
output_path = filepath where output figures should be saved.

USAGE: python volcano_plot_fragpipe.py </path/input_filename.txt> </path/output_foldername/>
'''

import os
import sys
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot

def volcano_plot(title, log2_fold_change, transformed_pvalues, annotations, SPO, present_only_in, output_path, x_axis_title="log2 fold change", y_axis_title="-log10 pvalue", point_radius=4):

    fig = go.Figure()
    fig.update_layout(
        title=title,
        xaxis_title=x_axis_title,
        yaxis_title=y_axis_title,
        paper_bgcolor='white',
        plot_bgcolor='white',
        shapes=[
            # Horizontal line for -log10(p-value) = 2
            dict(
                type='line',
                xref='paper', x0=0, x1=1,  # Line spans the entire x-axis
                yref='y', y0=2, y1=2,      # Line is fixed at y = 2
                line=dict(
                    color="grey",
                    width=1,
                    dash="dot",  # Dotted line style
                ),
            ),
            # Vertical line for log2 fold change = 2
            dict(
                type='line',
                xref='x', x0=2, x1=2,  # Line is fixed at x = 2
                yref='paper', y0=0, y1=1,  # Line spans the entire y-axis
                line=dict(
                    color="grey",
                    width=1,
                    dash="dot",  # Dotted line style
                ),
            ),
            # Vertical line for log2 fold change = -2
            dict(
                type='line',
                xref='x', x0=-2, x1=-2,  # Line is fixed at x = -2
                yref='paper', y0=0, y1=1,    # Line spans the entire y-axis
                line=dict(
                    color="grey",
                    width=1,
                    dash="dot",  # Dotted line style
                ),
            )
        ]
    )

    colors = []

    for i in range(len(log2_fold_change)):
        if "Only_" in str(present_only_in[i]):
            colors.append('#006400')  # dark green
        elif transformed_pvalues[i] > 2:
            if log2_fold_change[i] > 2:
                colors.append('#FFA500')  # yellow
            elif log2_fold_change[i] < -2:
                colors.append('#3f65d4')  # blue
            else:
                colors.append('rgba(150,150,150,0.5)')  # grey
        else:
            colors.append('rgba(150,150,150,0.5)')  # grey
            
    hover_texts = [f"SPO: {s} | Annotation: {a}" for s, a in zip(SPO, annotations)]

    fig.add_trace(
        go.Scattergl(
            x=log2_fold_change,
            y=transformed_pvalues,
            mode='markers',
            text=hover_texts,
            hovertemplate='%{text}<br>log2FC: %{x}<br>-log10(p): %{y}<extra></extra>',
            marker={
                'color': colors,
                'size': point_radius,
            }
        )
    )
    
    # Save the plot as an HTML file
    output_html = os.path.join(output_path, f"volcano_{title}.html")
    plot(fig, filename=output_html)

    # Save the plot as a PDF file
    output_pdf = os.path.join(output_path, f"volcano_{title}.pdf")
    fig.write_image(output_pdf)

if __name__ == "__main__": 
    if len(sys.argv) != 3:  # Check if there are exactly 2 arguments
        print("Usage: python volcano_plot.py </path/input_filename.txt> </path/output_foldername/>")
        sys.exit(1)
    
    inputfile = sys.argv[1]
    output_folder = sys.argv[2]

    # Load the input file into a DataFrame
    df = pd.read_csv(inputfile, sep='\t')

    # Filter out log2 fold change values > 15 or < -15
    filtered_df = df[(df.iloc[:, 14] <= 15) & (df.iloc[:, 14] >= -15)].reset_index(drop=True)

    # Extract relevant columns from the filtered DataFrame
    log2_fold_change = filtered_df.iloc[:, 14]
    transformed_pvalues = filtered_df.iloc[:, 17]
    SPO = filtered_df.iloc[:, 0]
    annotations = filtered_df.iloc[:, 1] 
    present_only_in = filtered_df["Present_Only_In"]

    title = os.path.basename(inputfile).replace(".txt", "")

    volcano_plot(
        title,
        log2_fold_change=log2_fold_change,
        transformed_pvalues=transformed_pvalues,
        annotations=annotations,
        SPO=SPO,
        present_only_in=present_only_in,
        output_path=output_folder
    )