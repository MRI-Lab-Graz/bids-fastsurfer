#!/usr/bin/env python3
import pandas as pd
import os

def augment_qdec(input_file, output_file, aseg_file=None):
    print(f"Reading {input_file}...")
    df = pd.read_csv(input_file, sep='\t')
    
    # If aseg file is provided, merge eTIV
    if aseg_file and os.path.exists(aseg_file):
        print(f"Reading eTIV from {aseg_file}...")
        # Read aseg, handling potential whitespace or tab delimiters
        try:
            aseg = pd.read_csv(aseg_file, sep='\t')
            if 'Measure:volume' not in aseg.columns: # Try whitespace if tab fails or looks wrong
                 aseg = pd.read_csv(aseg_file, delim_whitespace=True)
        except:
            aseg = pd.read_csv(aseg_file, delim_whitespace=True)
            
        # Standardize ID column
        id_col = 'Measure:volume' if 'Measure:volume' in aseg.columns else aseg.columns[0]
        
        if 'EstimatedTotalIntraCranialVol' in aseg.columns:
            etiv_data = aseg[[id_col, 'EstimatedTotalIntraCranialVol']].copy()
            etiv_data.columns = ['fsid', 'eTIV']
            
            # Fix ID mismatch: remove .long.template suffix if present
            # e.g., sub-01_ses-1.long.sub-01 -> sub-01_ses-1
            etiv_data['fsid'] = etiv_data['fsid'].apply(lambda x: x.split('.long.')[0] if '.long.' in str(x) else x)
            
            # Merge into QDEC
            df = pd.merge(df, etiv_data, on='fsid', how='left')
            print("Added column: eTIV")
        else:
            print("Warning: EstimatedTotalIntraCranialVol not found in aseg file")

    # Ensure group_beh_factor is treated as integer/string for mapping, not float
    # It might be read as int if no NaNs
    
    # Strategy 1: Duration (2w vs 4w vs Control)
    # 1 (alone_2w), 3 (group_2w) -> 2weeks
    # 2 (alone_4w), 4 (group_4w) -> 4weeks
    # 5 (control) -> control
    
    def map_duration(g):
        if g in [1, 3]: return "2weeks"
        if g in [2, 4]: return "4weeks"
        if g == 5: return "control"
        return "unknown"

    df['group_duration'] = df['group_beh_factor'].apply(map_duration)
    
    # Strategy 2: Context (Alone vs Group vs Control)
    # 1 (alone_2w), 2 (alone_4w) -> alone
    # 3 (group_2w), 4 (group_4w) -> group
    # 5 (control) -> control
    
    def map_context(g):
        if g in [1, 2]: return "alone"
        if g in [3, 4]: return "group"
        if g == 5: return "control"
        return "unknown"

    df['group_context'] = df['group_beh_factor'].apply(map_context)
    
    # Strategy 3: Intervention vs Control (Binary)
    # 1, 2, 3, 4 -> intervention
    # 5 -> control
    
    def map_binary(g):
        if g in [1, 2, 3, 4]: return "intervention"
        if g == 5: return "control"
        return "unknown"
        
    df['group_binary'] = df['group_beh_factor'].apply(map_binary)

    print("Added columns: group_duration, group_context, group_binary")
    print(df[['group_beh_factor', 'group_duration', 'group_context', 'group_binary']].drop_duplicates().sort_values('group_beh_factor'))
    
    print(f"Saving to {output_file}...")
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    augment_qdec("qdec.table.dat", "qdec.table.expanded.dat", "results/aseg.long.table")
