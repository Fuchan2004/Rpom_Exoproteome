#!/usr/bin/env python3

'''
This script extracts significant values (log2 fold change > 2 & -log10(p-value) > 2) from proteome data.

USAGE: python significant.py </path/input_filename> </path/output_folder/>
'''

import os
import sys
import pandas as pd
def significant_list(title, log2_fold_change, transformed_pvalues, annotations, accessions, condition, output_path):
    # Lists for standard over/underexpressed
    sign_plus, accession_sign_plus, fc_plus, pval_plus = [], [], [], []
    sign_minus, accession_sign_minus, fc_minus, pval_minus = [], [], [], []

    # Lists for extreme FC (zero in one condition)
    sign_plus_zero, accession_sign_plus_zero = [], []
    sign_minus_zero, accession_sign_minus_zero = [], []

    for i in range(len(log2_fold_change)):
        log2_fc = log2_fold_change[i]
        pval = transformed_pvalues[i]
        ann = annotations[i]
        acc = accessions[i]
        cond = condition[i]
        

        if pval > 2:
            if "Only_" in cond and log2_fc > 0:
                # Extremely overexpressed
                sign_plus_zero.append(ann)
                accession_sign_plus_zero.append(acc)
            if "Only_" in cond and log2_fc < 0:
                # Extremely underexpressed
                sign_minus_zero.append(ann)
                accession_sign_minus_zero.append(acc)
            if log2_fc > 2:
                # Standard overexpressed
                sign_plus.append(ann)
                accession_sign_plus.append(acc)
                fc_plus.append(log2_fc)
                pval_plus.append(pval)

                # Pad underexpressed side
                sign_minus.append(None)
                accession_sign_minus.append(None)
                fc_minus.append(None)
                pval_minus.append(None)

            elif log2_fc < -2:
                # Standard underexpressed
                sign_minus.append(ann)
                accession_sign_minus.append(acc)
                fc_minus.append(log2_fc)
                pval_minus.append(pval)

                # Pad overexpressed side
                sign_plus.append(None)
                accession_sign_plus.append(None)
                fc_plus.append(None)
                pval_plus.append(None)
            else:
                # Not significant FC
                sign_plus.append(None)
                accession_sign_plus.append(None)
                fc_plus.append(None)
                pval_plus.append(None)

                sign_minus.append(None)
                accession_sign_minus.append(None)
                fc_minus.append(None)
                pval_minus.append(None)
        else:
            # Not significant p-value
            sign_plus.append(None)
            accession_sign_plus.append(None)
            fc_plus.append(None)
            pval_plus.append(None)

            sign_minus.append(None)
            accession_sign_minus.append(None)
            fc_minus.append(None)
            pval_minus.append(None)

    # Create base DataFrame
    significant_df = pd.DataFrame({
        'Significant_Overexpressed': sign_plus,
        'Accession Number +': accession_sign_plus,
        'log2FC +': fc_plus,
        '-log10_adj_p +': pval_plus,
        'Significant_Underexpressed': sign_minus,
        'Accession Number -': accession_sign_minus,
        'log2FC -': fc_minus,
        '-log10_adj_p -': pval_minus
    })

    # Create extreme fold change dataframe
    zero_df = pd.DataFrame({
        'Zero_Condition_Overexpressed': sign_plus_zero + [None] * (max(len(sign_plus_zero), len(sign_minus_zero)) - len(sign_plus_zero)),
        'Accession Zero +': accession_sign_plus_zero + [None] * (max(len(sign_plus_zero), len(sign_minus_zero)) - len(accession_sign_plus_zero)),
        'Zero_Condition_Underexpressed': sign_minus_zero + [None] * (max(len(sign_plus_zero), len(sign_minus_zero)) - len(sign_minus_zero)),
        'Accession Zero -': accession_sign_minus_zero + [None] * (max(len(sign_plus_zero), len(sign_minus_zero)) - len(accession_sign_minus_zero))
    })

    # Combine
    final_df = pd.concat([significant_df, zero_df], axis=1)

    # Drop rows where all main values are None
    final_df.dropna(how='all', inplace=True)

    # Write to file
    output_filename = os.path.join(output_path, f"{title}_significant.txt")
    final_df.to_csv(output_filename, sep='\t', index=False)

    return final_df

if __name__ == "__main__": 
    if len(sys.argv) != 3:
        print("Usage: python significant.py </path/input_filename> </path/output_folder/>")
        sys.exit(1)
    
    inputfile = sys.argv[1]
    output_path = sys.argv[2]

    df = pd.read_csv(inputfile, sep='\t')

    log2_fold_change = df.iloc[:, 14]
    transformed_pvalues = df.iloc[:, 17]
    annotations = df.iloc[:, 1] 
    accessions = df.iloc[:, 0]
    condition = df.iloc[:, 9]

    title = os.path.basename(inputfile).replace(".txt", "")
    
    significant_list(
        title,
        log2_fold_change=log2_fold_change,
        transformed_pvalues=transformed_pvalues,
        annotations=annotations,
        accessions=accessions,
        condition=condition,
        output_path=output_path
    )