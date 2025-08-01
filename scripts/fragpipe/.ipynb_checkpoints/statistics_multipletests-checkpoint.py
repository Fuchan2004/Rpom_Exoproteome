#!/usr/bin/env python3

'''
This script generates statistics for FragPipe output data:
    1. Combines dataframes by Accession Number
    2. Calculates means and standard deviations
    3. Determines log2 fold changes
    4. Calculates p-values using user-specified statistical tests
    5. Exports the results

USAGE: python statistics.py <input_filename_1> <input_filename_2> <test_name> <multiple_test_correction>
    * Available tests: t-test, mannwhitneyu, wilcoxon, kruskal, f_oneway, ttest_rel, (default is t-test)
    * Multiple test correction: bonferroni, fdr_bh, None (default is None)

Explanation of tests: 
    A) Parametric Tests (Assume Normal Distribution)
        * t-test: Student t-test, Compares two independent samples.
        * ttest_rel: Paired t-test, Compares two related samples (e.g., before and after measurements).
        * f_oneway: ANOVA, Compares more than two groups for equality of means.
        * ttest_ind(..., equal_var=False): Welch's t-test, Variation of t-test, assumes unequal variances.

    B) Non-Parametric Tests (Assume no Normal Distribution)
        * mannwhitneyu: Mann-Whitney-U test, Compares two independent samples; detects differences in distribution.
        * wilcoxon: Wilcoxon signed rank, Compares two related samples; detects differences in medians.
        * kruskal: Kruskal Wallis test, Compares more than two groups (non-parametric version of ANOVA).
        
Recommendations for Multiple Testing Correction:
When performing many statistical tests simultaneously (e.g., calculating p-values for thousands of proteins in your proteomics data), the probability of obtaining false positives increases. This is known as the multiple testing problem. Multiple testing correction reduces this issue by adjusting the p-values to account for the number of tests performed, ensuring that the false discovery rate (FDR) or family-wise error rate (FWER) is controlled.

Common Multiple Testing Correction Methods:
    * bonferroni: Bonferroni correction (most stringent), Controls the FWER but can be overly conservative.
    * fdr_bh: Benjamini-Hochberg correction, Correction most widely used in omics data, Adjusts for the False Discovery Rate (FDR). More powerful than Bonferroni for large datasets.
    * fdr_by: Benjamini-Yekutieli correction, A more conservative version of fdr_bh, applicable when test results are dependent.
    * holm: Holm correction, Less conservative than Bonferroni.
    * sidak: correction: Similar to Bonferroni but slightly less stringent.
'''

import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu, wilcoxon, kruskal, f_oneway
from statsmodels.stats.multitest import multipletests
import math

# Map test names to functions
def get_test_function(test_name):
    test_mapping = {
        't-test': ttest_ind,
        'welch': 'welch',
        'mannwhitneyu': mannwhitneyu,
        'wilcoxon': wilcoxon,
        'kruskal': kruskal,
        'f_oneway': f_oneway
    }
    return test_mapping.get(test_name.lower(), ttest_ind)  # Default to t-test


def merge_and_format(df1, df2, suf_1='', suf_2=''):
    """ Merges and formats the two input dataframes. """
    merged_df = pd.merge(df1, df2, on="Accession Number", how="outer", suffixes=(suf_1, suf_2))

    # Replace 0 in annotations with "Unknown"
    if f'Annotation{suf_1}' in merged_df.columns:
        merged_df[f'Annotation{suf_1}'] = merged_df[f'Annotation{suf_1}'].replace(0, "Unknown")
    if f'Annotation{suf_2}' in merged_df.columns:
        merged_df[f'Annotation{suf_2}'] = merged_df[f'Annotation{suf_2}'].replace(0, "Unknown")

    return merged_df


def calculate_log2_fold_change(df, suf_1, suf_2):
    """ Calculate log2 fold change """
    # Check if any value in the Row_Average is 0, and replace it with a small number
    if (df[f'Row_Average_{suf_2}'] == 0).any():  # Check if any value is 0
        df[f'Row_Average_{suf_2}'] = df[f'Row_Average_{suf_2}'].replace(0, 10E-15)
    if (df[f'Row_Average_{suf_1}'] == 0).any():  # Check if any value is 0
        df[f'Row_Average_{suf_1}'] = df[f'Row_Average_{suf_1}'].replace(0, 10E-15)
    log2_fc = np.log2((df[f'Row_Average_{suf_1}']) / (df[f'Row_Average_{suf_2}'])).replace([np.inf, -np.inf], np.nan).dropna()
    return log2_fc


def calculate_pvalues(df, suf_1, suf_2, test_func):
    """ Perform row-wise tests and return raw p-values """
    p_values = []

    for _, row in df.iterrows():
        try:
            group1 = row[[f'1_{suf_1}', f'2_{suf_1}', f'3_{suf_1}']].astype(float).values
            group2 = row[[f'1_{suf_2}', f'2_{suf_2}', f'3_{suf_2}']].astype(float).values

            # Replace with small values if either group contains only zeros
            if np.all(group1 == 0): 
                group1 = [10E-15, 10E-15, 10E-15,]
            elif np.all(group2 == 0):
                group2 = [10E-15, 10E-15, 10E-15,]  
                #   p_values.append(1.0)
            #else:
                # Handle Welch test separately
            if test_func == 'welch':
                _, p_val = ttest_ind(group1, group2, equal_var=False, nan_policy='omit')
            elif test_func in {ttest_ind, mannwhitneyu, wilcoxon, kruskal}:
                _, p_val = test_func(group1, group2, nan_policy='omit')
            elif test_func == f_oneway:
                p_val = test_func(group1, group2).pvalue
            else:
                p_val = 1.0  # Default to 1 if invalid test

            p_values.append(p_val)

        except ValueError:
            p_values.append(1.0)

    return np.array(p_values)


def multiple_testing_correction(p_values, method='None'):
    """ Applies multiple testing correction and returns both corrected p-values and their negative log-transformed values. """
    
    # Valid methods for multiple testing correction
    valid_methods = ['bonferroni', 'holm', 'sidak', 'fdr_bh', 'fdr_by']
    
    # If no correction method is chosen, return raw p-values
    if method == "None":
        print("No multiple testing correction applied.")
        return p_values, -np.log10(p_values)  # Return raw p-values and their log-transformed values

    # If an invalid method is provided, print a warning and use no correction
    if method not in valid_methods:
        print(f"Invalid correction method '{method}', using 'None' by default.")
        method = 'None'
        return p_values, -np.log10(p_values)

    # Apply multiple testing correction
    _, corrected_pvals, _, _ = multipletests(p_values, method=method)

    # Return both corrected p-values and their negative log-transformed values
    return corrected_pvals, -np.log10(corrected_pvals)

def valid_triplicate(values):
    nan_count = sum(math.isnan(v) for v in values)
    return nan_count == 0 or nan_count == 3

def is_all_zero(values):
    return all(v == 0 for v in values)

def valid_row(row, suf_1, suf_2):
    group1 = [row[c] for c in [f'1_{suf_1}', f'2_{suf_1}', f'3_{suf_1}']]
    group2 = [row[c] for c in [f'1_{suf_2}', f'2_{suf_2}', f'3_{suf_2}']]

    # drop row if both groups are all zero
    if is_all_zero(group1) and is_all_zero(group2):
        return False

    return valid_triplicate(group1) and valid_triplicate(group2)


def process_files(file_1, file_2, test_name, correction_method='None'):
    """ Main function to process and merge files, calculate stats, and export the result. """
    
    test_func = get_test_function(test_name)

    # Extract metadata from filenames
    filename_1 = os.path.basename(file_1)
    filename_2 = os.path.basename(file_2)

    strain = filename_1.split('_')[1]
    medium_1 = filename_1.split('_')[2]
    medium_2 = filename_2.split('_')[2]
    phase_1 = filename_1.split('_')[3]
    phase_2 = filename_2.split('_')[3]

    suf_1 = phase_1 if medium_1 == medium_2 else medium_1
    suf_2 = phase_2 if medium_1 == medium_2 else medium_2

    # Read files into pandas dataframes
    df1 = pd.read_csv(file_1, sep='\t')
    df2 = pd.read_csv(file_2, sep='\t')

    # Convert measurements to numeric, coercing errors
    df1.iloc[:, 3:6] = df1.iloc[:, 3:6].apply(pd.to_numeric, errors='coerce')
    df2.iloc[:, 3:6] = df2.iloc[:, 3:6].apply(pd.to_numeric, errors='coerce')

    # Rename replicate columns
    for i in range(1, 4):
        df1.rename(columns={df1.columns[i + 1]: f'{i}_{suf_1}'}, inplace=True)
        df2.rename(columns={df2.columns[i + 1]: f'{i}_{suf_2}'}, inplace=True)

    # Merge dataframes
    combined_df = merge_and_format(df1, df2, suf_1, suf_2)

    # Drop rows with missing values in any of the replicate columns
    print("Before filtering:", len(combined_df))
    combined_df = combined_df[combined_df.apply(lambda row: valid_row(row, suf_1, suf_2), axis=1)]
    print("After filtering:", len(combined_df))
    combined_df.fillna(0, inplace=True)

    def presence_annotation(row, suf_1, suf_2):
        group1 = [row[f'1_{suf_1}'], row[f'2_{suf_1}'], row[f'3_{suf_1}']]
        group2 = [row[f'1_{suf_2}'], row[f'2_{suf_2}'], row[f'3_{suf_2}']]
    
        if is_all_zero(group1) and not is_all_zero(group2):
            return f'Only_{suf_2}'
        elif is_all_zero(group2) and not is_all_zero(group1):
            return f'Only_{suf_1}'
        else:
            return 'Both'

    combined_df['Present_Only_In'] = combined_df.apply(lambda row: presence_annotation(row, suf_1, suf_2), axis=1)


    # Calculate row-wise averages and standard deviations
    combined_df[f'Row_Average_{suf_1}'] = combined_df[[f'1_{suf_1}', f'2_{suf_1}', f'3_{suf_1}']].mean(axis=1)
    combined_df[f'STD_{suf_1}'] = combined_df[[f'1_{suf_1}', f'2_{suf_1}', f'3_{suf_1}']].std(axis=1)

    combined_df[f'Row_Average_{suf_2}'] = combined_df[[f'1_{suf_2}', f'2_{suf_2}', f'3_{suf_2}']].mean(axis=1)
    combined_df[f'STD_{suf_2}'] = combined_df[[f'1_{suf_2}', f'2_{suf_2}', f'3_{suf_2}']].std(axis=1)

    # Calculate log2 fold change and p-values
    combined_df['Log2_Fold_Change'] = calculate_log2_fold_change(combined_df, suf_1, suf_2)

    raw_pvals = calculate_pvalues(combined_df, suf_1, suf_2, test_func)

    if correction_method != 'None':
        corrected_pvals, neglog_pvals = multiple_testing_correction(raw_pvals, correction_method)
        combined_df['Adj_P_Values'] = corrected_pvals
        combined_df['-log10(Adj_P_Values)'] = neglog_pvals
    else:
        combined_df['Adj_P_Values'] = raw_pvals
        combined_df['-log10(Adj_P_Values)'] = -np.log10(raw_pvals)


    # Output directory and filename formatting
    input_directory = os.path.dirname(file_1)
    output_directory =  'output/proteomics/'
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Calculate log10 transformed p-values for raw and adjusted p-values
    if correction_method != 'None':
        combined_df['Log10_Adj_P_Values'] = -np.log10(combined_df['Adj_P_Values'])
    else:
        combined_df['Log10_P_Values'] = -np.log10(combined_df['P_Values'])

    if medium_1 == medium_2:
        output_filename = os.path.join(output_directory, f'{strain}_{medium_1}_{phase_1}VS{phase_2}_{test_name}_{correction}.txt')
    else:
        output_filename = os.path.join(output_directory, f'{strain}_{phase_1}_{medium_1}VS{medium_2}_{test_name}_{correction}.txt')

    combined_df.to_csv(output_filename, sep='\t', index=False)
    print(f"Data saved to {output_filename}")


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        print("USAGE: python statistics.py <input_file1> <input_file2> <test_name> <correction_method>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    test = sys.argv[3] if len(sys.argv) > 3 else 't-test'
    correction = sys.argv[4] if len(sys.argv) > 4 else 'None'

    process_files(file1, file2, test, correction)

    # Handle 3 or 4 arguments
    if len(sys.argv) == 3:
        # Only two files → t-test with no correction
        pass
    elif len(sys.argv) == 4:
        arg = sys.argv[3].lower()
    
        # Check if the third argument is a valid correction method
        valid_corrections = {'bonferroni', 'fdr_bh', 'fdr_by', 'holm', 'sidak'}
    
        if arg in valid_corrections:
            correction = arg  # Apply correction with t-test
        else:
            test_name = arg  # Use specified test without correction
    elif len(sys.argv) == 5:
        # Both test and correction are specified
        test_name = sys.argv[3].lower()
        correction = sys.argv[4].lower()

    # Display the run details
    print(f"Statistics run with {test_name} and {correction}")

    # Run the analysis
    process_files(file1, file2, test_name, correction)