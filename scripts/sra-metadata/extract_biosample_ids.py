#!/usr/bin/env python3
"""
Extract BioSample IDs from NCBI SRA_Accessions.tab using column names.

This script streams the large SRA_Accessions.tab file (29 GB) from NCBI FTP,
parses the header to find column positions dynamically, and extracts only the specified
columns for matching SRA IDs. This approach is robust to future changes in column ordering.

Usage:
    # Download and stream from NCBI FTP (recommended - no local storage of 29 GB file)
    python extract_biosample_ids.py \\
        --sra-ids input_sra_ids.txt \\
        --output sra_biosample_map.txt \\
        --ncbi-url ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab \\
        --columns Accession BioSample

Educational notes:
    - Streaming approach: Processes 29 GB file without storing it locally
    - Header-based: Uses column names (not positions) for robustness
    - Set-based lookup: Loads SRA IDs into set for O(1) lookup time
    - Progress reporting: Shows matched IDs every 10,000 rows

Author: SRA metadata extraction pipeline
"""

import sys
import argparse
import subprocess
from typing import Set, List, Iterator

def load_sra_ids(filename: str) -> Set[str]:
    """
    Load SRA IDs from file into a set for fast O(1) lookup.

    Args:
        filename: Path to file with SRA IDs (one per line)

    Returns:
        Set of SRA ID strings

    Educational note:
        Using a set instead of a list provides O(1) lookup time vs O(n),
        which is critical when processing millions of rows from SRA_Accessions.tab
    """
    with open(filename, 'r') as f:
        return {line.strip() for line in f if line.strip()}

def stream_ncbi_file(ncbi_url: str) -> Iterator[str]:
    """
    Stream lines from NCBI FTP using wget.

    Args:
        ncbi_url: URL to NCBI SRA_Accessions.tab file

    Yields:
        Lines from the file

    Educational note:
        Uses subprocess to stream wget output line-by-line without
        storing the entire 29 GB file in memory or on disk.
    """
    process = subprocess.Popen(
        ['wget', '-qO-', ncbi_url],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1  # Line-buffered
    )

    for line in process.stdout:
        yield line

    # Check for errors
    process.wait()
    if process.returncode != 0:
        stderr = process.stderr.read()
        raise RuntimeError(f"wget failed with return code {process.returncode}: {stderr}")

def extract_columns_by_name(sra_ids: Set[str],
                            output_file: str,
                            column_names: List[str],
                            input_stream: Iterator[str]):
    """
    Stream input, extract specified columns by name for matching SRA IDs.

    This function:
    1. Reads the header line to find column positions
    2. Validates all requested columns exist
    3. Streams data rows, filtering to matching SRA IDs
    4. Writes only requested columns to output

    Args:
        sra_ids: Set of SRA IDs to extract
        output_file: Path to output TSV file
        column_names: List of column names to extract (must include "Accession")
        input_stream: Iterator yielding lines from the input file

    Educational notes:
        - Column indices are found dynamically from header (robust to reordering)
        - Streaming processing: Never loads entire 29 GB file into memory
        - Progress reporting: Updates every 10,000 matches
    """
    # Step 1: Read header from input stream and find column positions
    header_line = next(input_stream)
    header_cols = header_line.strip().split('\t')

    # Step 2: Validate requested columns exist and get their indices
    try:
        col_indices = {name: header_cols.index(name) for name in column_names}
    except ValueError as e:
        # If a requested column doesn't exist, show error and available columns
        print(f"Error: Column not found in header: {e}", file=sys.stderr)
        print(f"Available columns: {', '.join(header_cols)}", file=sys.stderr)
        sys.exit(1)

    # Step 3: Validate that "Accession" column is requested (needed for filtering)
    if 'Accession' not in col_indices:
        print("Error: 'Accession' must be in column_names for filtering",
              file=sys.stderr)
        sys.exit(1)

    # Step 4: Write header to output file and process data
    with open(output_file, 'w') as out_f:
        # Write output header
        out_f.write('\t'.join(column_names) + '\n')

        # Initialize counters for progress reporting
        matched = 0      # Number of SRA IDs found
        processed = 0    # Total rows processed

        # Step 5: Stream and process data lines
        for line in input_stream:
            processed += 1
            fields = line.strip().split('\t')

            # Get the Accession (SRA ID) for this row
            accession = fields[col_indices['Accession']]

            # Step 6: Filter - only process if SRA ID is in our input set
            if accession in sra_ids:
                # Extract requested columns in the order specified
                row_data = [fields[col_indices[name]] for name in column_names]
                out_f.write('\t'.join(row_data) + '\n')
                matched += 1

                # Progress reporting every 10,000 matches
                if matched % 10000 == 0:
                    print(f"Matched {matched:,} SRA IDs (processed {processed:,} rows)...",
                          file=sys.stderr)

    # Final summary
    print(f"\nComplete!", file=sys.stderr)
    print(f"  Matched: {matched:,} out of {len(sra_ids):,} requested SRA IDs",
          file=sys.stderr)
    print(f"  Processed: {processed:,} total rows from SRA_Accessions.tab",
          file=sys.stderr)
    print(f"  Output: {output_file}",
          file=sys.stderr)

    # Warn if not all SRA IDs were found
    if matched < len(sra_ids):
        print(f"\n  Warning: {len(sra_ids) - matched:,} SRA IDs not found in SRA_Accessions.tab",
              file=sys.stderr)

def main():
    """Parse arguments and run column extraction."""
    parser = argparse.ArgumentParser(
        description='Extract columns from SRA_Accessions.tab by column name',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download and stream from NCBI FTP (recommended)
  python extract_biosample_ids.py \\
      --sra-ids my_sra_ids.txt \\
      --output sra_biosample_map.txt \\
      --ncbi-url ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab \\
      --columns Accession BioSample

Column names available in SRA_Accessions.tab (as of 2024):
  Accession, Submission, Status, Updated, Published, Received, Type,
  Center, Visibility, Alias, Experiment, Sample, Study, Loaded, Spots,
  Bases, Md5sum, BioSample, BioProject, ReplacedBy
        """)

    parser.add_argument('--sra-ids', required=True,
                       help='File with SRA IDs to extract (one per line)')
    parser.add_argument('--output', required=True,
                       help='Output TSV file with extracted columns')
    parser.add_argument('--ncbi-url', required=True,
                       help='URL to NCBI SRA_Accessions.tab file')
    parser.add_argument('--columns', nargs='+', required=True,
                       help='Column names to extract (must include "Accession")')

    args = parser.parse_args()

    # Load SRA IDs into memory as a set (for fast lookup)
    print(f"Loading SRA IDs from {args.sra_ids}...", file=sys.stderr)
    sra_ids = load_sra_ids(args.sra_ids)
    print(f"Loaded {len(sra_ids):,} SRA IDs to match", file=sys.stderr)

    # Stream from NCBI and extract columns
    print(f"\nExtracting columns: {', '.join(args.columns)}", file=sys.stderr)
    print(f"Streaming from {args.ncbi_url}...", file=sys.stderr)

    input_stream = stream_ncbi_file(args.ncbi_url)
    extract_columns_by_name(sra_ids, args.output, args.columns, input_stream)

if __name__ == "__main__":
    main()
