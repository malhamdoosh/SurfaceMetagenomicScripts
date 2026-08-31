import os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def main():
    # 1. Setup Command Line Arguments
    parser = argparse.ArgumentParser(description="Generate clinical validation heatmap from Kraken and Assay data.")
    
    # Input files
    parser.add_argument('--targets', type=str, default='target_list.txt', help='Path to target list TSV')
    parser.add_argument('--assay', type=str, default='assay_results.txt', help='Path to assay results TSV')
    parser.add_argument('--reports_dir', type=str, default='.', help='Directory containing _report.txt files')
    
    # Output names
    parser.add_argument('--out_matrix', type=str, default='pathogen_presence_absence_matrix.csv', help='Output filename for PA matrix')
    parser.add_argument('--out_png', type=str, default='heatmap.png', help='Output filename for the image')

    args = parser.parse_args()

    # 2. Load Target List
    if not os.path.exists(args.targets):
        print(f"Error: {args.targets} not found.")
        return

    targets_df = pd.read_csv(
        args.targets, 
        sep='\t', 
        names=['rank', 'TaxID', 'name_kraken_report', 'target_name_assay'],
        dtype={'TaxID': str},
        skiprows=1
    )
    targets_df['TaxID'] = targets_df['TaxID'].str.strip()

    # 3. Process Kraken Reports
    report_files = sorted([f for f in os.listdir(args.reports_dir) if f.endswith('_report.txt')])
    if not report_files:
        print(f"No _report.txt files found in {args.reports_dir}")
        return

    sample_columns = []
    for file in report_files:
        sample_name = os.path.splitext(file)[0].split("_")[0]
        file_path = os.path.join(args.reports_dir, file)
        try:
            report_data = pd.read_csv(
                file_path, sep='\t', header=None, usecols=[1, 4], 
                names=['Reads', 'TaxID'], dtype={'TaxID': str}
            )
            report_data['TaxID'] = report_data['TaxID'].str.strip()
            merged = pd.merge(targets_df[['TaxID']], report_data, on='TaxID', how='left')
            sample_col = merged[['Reads']].rename(columns={'Reads': sample_name})
            sample_columns.append(sample_col)
        except Exception as e:
            print(f"Error processing {file}: {e}")

    # Save Grand Presence/Absence Matrix (Required for raw data verification)
    pa_matrix_numeric = pd.concat([targets_df] + sample_columns, axis=1)
    pa_matrix_numeric.fillna("NA").to_csv(args.out_matrix, index=False)

    # 4. Process Assay Results
    if not os.path.exists(args.assay):
        print(f"Error: {args.assay} not found.")
        return

    assay_results = pd.read_csv(args.assay, sep='\t', names=['sample_name', 'assay_targets'])
    all_potential_targets = sorted(targets_df['target_name_assay'].unique())

    # 5. Generate Heatmap Logic (Internal DataFrame only)
    heatmap_rows = []
    for _, row in assay_results.iterrows():
        s_name = str(row['sample_name'])
        assay_positives = [t.strip() for t in str(row['assay_targets']).split(',')]
        row_data = {'Sample_ID': s_name}
        
        for target in all_potential_targets:
            is_assay_pos = target in assay_positives
            if s_name in pa_matrix_numeric.columns:
                target_slice = pa_matrix_numeric[pa_matrix_numeric['target_name_assay'] == target][s_name]
                
                if target_slice.isna().all():
                    res = "NA"
                else:
                    reads = pd.to_numeric(target_slice, errors='coerce').fillna(0).sum()
                    has_reads = reads >= 1
                    if is_assay_pos:
                        res = 'TP' if has_reads else 'FN'
                    else:
                        res = 'FP' if has_reads else 'TN'
            else:
                res = "NA"
            row_data[target] = res
        heatmap_rows.append(row_data)

    heatmap_df = pd.DataFrame(heatmap_rows)

    # 6. Visualization
    plot_df = heatmap_df.set_index('Sample_ID')
    map_dict = {'NA': 0, 'TN': 1, 'TP': 2, 'FP': 3, 'FN': 4}
    numeric_df = plot_df.replace(map_dict).apply(pd.to_numeric, errors='coerce').fillna(0)

    # Colors: Gray, Light Green, Dark Green, Orange, Red
    colors = ["#d1d1d1", "#a8e6cf", "#2e7d32", "#ffc107", "#f44336"]
    labels = ["NA (No Data)", "TN (True Negative)", "TP (True Positive)", "FP (Unexpected Hit)", "FN (Missing Hit)"]

    fig, ax = plt.subplots(figsize=(14, 8))
    sns.heatmap(numeric_df, cmap=colors, cbar=False, linewidths=.5, linecolor='white', ax=ax)

    # Legend setup
    patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(colors))]
    ax.legend(handles=patches, bbox_to_anchor=(1.02, 1), loc='upper left', title="Validation Key")
    
    ax.set_title("Metagenomics vs. QIAstat-DX assay", fontsize=16, pad=20)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Save PNG Only
    plt.savefig(args.out_png, dpi=300)

    print(f"\nProcessing complete!")
    print(f"Raw data saved to: {args.out_matrix}")
    print(f"Visual heatmap saved to: {args.out_png}")

if __name__ == "__main__":
    main()
