#!/usr/bin/env python3
import os
import sys
import csv
import argparse
import re

def parse_aseg_stats(stats_file):
    data = {
        'Left-Caudate': None,
        'Left-Putamen': None,
        'Left-Pallidum': None,
        'Right-Caudate': None,
        'Right-Putamen': None,
        'Right-Pallidum': None,
        'eTIV': None
    }
    
    try:
        with open(stats_file, 'r') as f:
            for line in f:
                line = line.strip()
                # Parse eTIV from header
                if line.startswith('# Measure EstimatedTotalIntraCranialVol') or line.startswith('# Measure eTIV'):
                    parts = line.split(',')
                    if len(parts) >= 4:
                        data['eTIV'] = parts[3].strip()
                
                # Parse table rows
                if not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 5:
                        struct_name = parts[4]
                        volume = parts[3]
                        if struct_name in data:
                            data[struct_name] = volume
                            
    except Exception as e:
        print(f"Error reading {stats_file}: {e}", file=sys.stderr)
        return None

    return data

def main():
    parser = argparse.ArgumentParser(description="Extract subcortical volumes (Caudate, Putamen, Pallidum) and eTIV from FreeSurfer/FastSurfer output.")
    parser.add_argument("subjects_dir", help="Path to the subjects directory (FreeSurfer/FastSurfer output folder)")
    parser.add_argument("output_file", help="Path to the output TSV file")
    parser.add_argument("--spss", help="Path to output SPSS-friendly CSV file", default=None)
    args = parser.parse_args()

    subjects_dir = args.subjects_dir
    output_file = args.output_file

    if not os.path.isdir(subjects_dir):
        print(f"Error: {subjects_dir} is not a directory.", file=sys.stderr)
        sys.exit(1)

    # Columns for the output file
    columns = ['subject_id_session', 'Left-Caudate', 'Left-Putamen', 'Left-Pallidum', 'Right-Caudate', 'Right-Putamen', 'Right-Pallidum', 'eTIV']
    spss_columns = ['Subject', 'Session'] + columns
    
    results = []

    # Iterate over subdirectories in subjects_dir
    # We assume any directory containing stats/aseg.stats is a subject/session directory
    print(f"Scanning {subjects_dir} for stats/aseg.stats files...")
    
    subjects = sorted([d for d in os.listdir(subjects_dir) if os.path.isdir(os.path.join(subjects_dir, d))])
    
    # Filter out cross-sectional runs if longitudinal run exists
    filtered_subjects = []
    for s in subjects:
        if ".long." not in s:
            # Check if a longitudinal version exists
            # We look for any folder that starts with s + ".long."
            has_long = any(other.startswith(s + ".long.") for other in subjects)
            if has_long:
                continue
        filtered_subjects.append(s)
    subjects = filtered_subjects

    for subject_id in subjects:
        stats_path = os.path.join(subjects_dir, subject_id, 'stats', 'aseg.stats')
        
        if os.path.isfile(stats_path):
            print(f"Processing {subject_id}...")
            stats_data = parse_aseg_stats(stats_path)
            
            if stats_data:
                # Extract Subject and Session
                subj_match = re.search(r'sub-([a-zA-Z0-9]+)_ses-([a-zA-Z0-9]+)', subject_id)
                subject_code = subj_match.group(1) if subj_match else subject_id
                session_code = subj_match.group(2) if subj_match else "NA"

                row = {
                    'subject_id_session': subject_id,
                    'Subject': subject_code,
                    'Session': session_code,
                    'Left-Caudate': stats_data['Left-Caudate'],
                    'Left-Putamen': stats_data['Left-Putamen'],
                    'Left-Pallidum': stats_data['Left-Pallidum'],
                    'Right-Caudate': stats_data['Right-Caudate'],
                    'Right-Putamen': stats_data['Right-Putamen'],
                    'Right-Pallidum': stats_data['Right-Pallidum'],
                    'eTIV': stats_data['eTIV']
                }
                results.append(row)
        else:
            # Check if it might be a BIDS structure or deeper nesting? 
            # For now, assume flat subject/session directories as is typical in FS output.
            pass

    if not results:
        print("No aseg.stats files found in immediate subdirectories.", file=sys.stderr)
        sys.exit(1)

    print(f"Writing results to {output_file}...")
    try:
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            writer.writerows(results)
        print("Done.")
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)

    if args.spss:
        print(f"Writing SPSS file to {args.spss}...")
        try:
            import pandas as pd
            import pyreadstat
            
            # Create DataFrame
            df = pd.DataFrame(results)
            
            # Ensure numeric columns are actually numeric
            numeric_cols = ['Left-Caudate', 'Left-Putamen', 'Left-Pallidum', 'Right-Caudate', 'Right-Putamen', 'Right-Pallidum', 'eTIV']
            for col in numeric_cols:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # Reorder columns to put Subject/Session first if they exist
            cols = list(df.columns)
            if 'Subject' in cols and 'Session' in cols:
                # Move Subject and Session to the front
                cols.insert(0, cols.pop(cols.index('Subject')))
                cols.insert(1, cols.pop(cols.index('Session')))
                df = df[cols]

            # Sanitize column names for SPSS (no hyphens)
            df.columns = [c.replace('-', '_') for c in df.columns]

            pyreadstat.write_sav(df, args.spss)
            print("Done writing SPSS file.")
            
        except ImportError:
            print("Error: pandas and pyreadstat are required for SPSS output. Please install them (pip install pandas pyreadstat).", file=sys.stderr)
        except Exception as e:
            print(f"Error writing SPSS file: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
