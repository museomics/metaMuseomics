import os
import csv
import sys
import pandas
import pathlib

def get_ids(file_path, column_name):
    # Reading in file 
    df = pandas.read_csv(file_path)

    # Extract IDS
    ids = df[column_name]

    return ids

def find_files(project_dir, ids):
    # Initialize a dictionary to store the results
    results = {}

    # Loop through each id
    for id in ids:
        dir = pathlib.Path(project_dir)
        files = [str(file) for file in dir.glob(f"**/*{id}*trimmed*.f*q.gz")]
        files = sorted(files)
       
        if len(files)==0:
            print(f"No files found for ID: {id}")  
            continue
        if len(files)==1:
            print(f"Your read files are unpaired for ID {id}")
            continue
        if len(files)>2:
            print(f"[WARNING] Multiple files found for ID '{id}' expected exactly 2. Files: {files}")
            continue
        
        r1_path = files[0]
        r2_path = files[1]
        # this gives relative file path for full path:  str(files[0].resolve())
        
        if r1_path and r2_path:
            # Process ID
            results[id] = (r1_path, r2_path)

    return results

def write_to_csv(file_path, column_name, results, output_filename):
    try:
        # Read the original input file
        df = pandas.read_csv(file_path)

        # Add empty columns for forward/reverse
        df["forward"] = ""
        df["reverse"] = ""

        # Fill in paths where IDs match
        for id, (r1_path, r2_path) in results.items():
            df.loc[df[column_name] == id, "forward"] = r1_path
            df.loc[df[column_name] == id, "reverse"] = r2_path

        # Write back out (keep same delimiter)
        df.to_csv(output_filename, index=False)
        print(f"[SUCCESS] Updated file written to: {output_filename}")

    except Exception as e:
        print(f"[ERROR] Failed to append R1/R2 columns to output CSV: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <project_dir> <sample_info_file> <column_name>")
        sys.exit(1)

    project_dir = sys.argv[1]
    output_filename = os.getcwd() + "_samples_out.csv"
    file_path = sys.argv[2]
    column_name = sys.argv[3]
    ids = get_ids(file_path, column_name)
    #print("IDs", ids)
    files_info = find_files(project_dir, ids)
    write_to_csv(file_path, column_name, files_info, output_filename)
    print(f"CSV file '{output_filename}' created successfully.")
