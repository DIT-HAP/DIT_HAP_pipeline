#!/usr/bin/env python3
"""
Download KEGG BRITE JSON files for S. pombe from the KEGG database.
Reads BRITE codes from pombe_kegg_brite.xlsx and downloads corresponding JSON files.

Usage:
    # Download all BRITE codes from Excel file
    python download_pombe_brite.py
    
    # Download specific BRITE code(s)
    python download_pombe_brite.py 03032
    python download_pombe_brite.py 03032 03011 01000
    
    # Specify organism (default: spo for S. pombe)
    python download_pombe_brite.py --organism spo 03032
"""

import pandas as pd
import requests
from pathlib import Path
import time
from typing import Optional, List
import json
import argparse
import sys


def download_brite_json(brite_code: str, output_dir: Path, organism: str = "spo") -> Optional[Path]:
    """
    Download BRITE JSON file from KEGG database.
    
    Args:
        brite_code: BRITE code (e.g., "03032")
        output_dir: Directory to save the JSON file
        organism: Organism code (default: "spo" for S. pombe)
    
    Returns:
        Path to downloaded file if successful, None otherwise
    """
    # Construct the URL
    url = f"https://www.kegg.jp/kegg-bin/download_htext?htext={organism}{brite_code}&format=json&filedir=kegg/brite/{organism}"
    
    # Create output filename
    output_file = output_dir / f"{organism}{brite_code}.json"
    
    try:
        print(f"Downloading {organism}{brite_code}...")
        print(f"  URL: {url}")
        response = requests.get(url, timeout=10)
        print(f"  Status: {response.status_code}")
        response.raise_for_status()
        
        # Check content type and length
        content_type = response.headers.get('content-type', '')
        content_length = len(response.content)
        print(f"  Content-Type: {content_type}, Length: {content_length} bytes")
        
        # Check if we got valid JSON
        try:
            data = response.json()
            # Check if it's an empty or error response
            if not data or (isinstance(data, dict) and not data.get('children')):
                print(f"  ⚠ Empty or invalid BRITE hierarchy for {brite_code}")
                return None
            
            # Save to file
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            print(f"  ✓ Saved to {output_file}")
            return output_file
        except json.JSONDecodeError as je:
            print(f"  ✗ Invalid JSON response for {brite_code}: {je}")
            # Save the raw response for debugging
            debug_file = output_dir / f"{organism}{brite_code}_debug.txt"
            with open(debug_file, 'w', encoding='utf-8') as f:
                f.write(response.text[:1000])  # First 1000 chars
            print(f"  Debug: First 1000 chars saved to {debug_file}")
            return None
            
    except requests.exceptions.Timeout:
        print(f"  ✗ Timeout downloading {brite_code} (>10s)")
        return None
    except requests.exceptions.RequestException as e:
        print(f"  ✗ Error downloading {brite_code}: {e}")
        return None


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Download KEGG BRITE JSON files for S. pombe',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Download all BRITE codes from Excel file
  %(prog)s
  
  # Download specific BRITE code(s)
  %(prog)s 03032
  %(prog)s 03032 03011 01000
  
  # Specify organism
  %(prog)s --organism sce 03032
        '''
    )
    
    parser.add_argument(
        'codes',
        nargs='*',
        help='BRITE code(s) to download (e.g., 03032, 03011). If not provided, downloads all codes from Excel file.'
    )
    
    parser.add_argument(
        '-o', '--organism',
        default='spo',
        help='Organism code (default: spo for S. pombe)'
    )
    
    parser.add_argument(
        '--output-dir',
        type=Path,
        help='Output directory for JSON files (default: ./brite_json)'
    )
    
    return parser.parse_args()


def get_brite_codes_from_excel(excel_file: Path) -> List[str]:
    """Load BRITE codes from Excel file."""
    print(f"Loading BRITE codes from {excel_file}...")
    df = pd.read_excel(excel_file, index_col=[0])
    
    # Ensure BRITE codes are 5 characters with leading zeros
    if code_column := next((col for col in df.columns if any(kw in col.lower() for kw in ['code', 'brite', 'id', 'number'])), df.columns[0] if len(df.columns) > 0 else None):
        df[code_column] = df[code_column].astype(int).astype(str).str.strip().str.zfill(5)
    
    # Display the table
    print(f"\nFound {len(df)} BRITE entries:")
    print(df.to_string(index=False))
    print()
    
    # Determine which column contains the BRITE codes
    code_column = None
    for col in df.columns:
        col_lower = col.lower()
        if any(keyword in col_lower for keyword in ['code', 'brite', 'id', 'number']):
            code_column = col
            break
    
    if code_column is None and len(df.columns) > 0:
        code_column = df.columns[0]
    
    print(f"Using column '{code_column}' for BRITE codes\n")
    
    # Extract BRITE codes
    brite_codes = []
    for idx, row in df.iterrows():
        brite_code = str(row[code_column]).strip()
        
        # Extract numeric code if needed (e.g., "ko03032" -> "03032")
        if brite_code.startswith('ko'):
            brite_code = brite_code[2:]
        elif brite_code.startswith('spo'):
            brite_code = brite_code[3:]
        
        brite_codes.append(brite_code)
    
    return brite_codes


def normalize_brite_code(code: str) -> str:
    """Normalize BRITE code by removing prefix and ensuring proper format."""
    code = code.strip()
    
    # Remove common prefixes
    if code.startswith('ko'):
        code = code[2:]
    elif code.startswith('spo'):
        code = code[3:]
    elif code.startswith('br:'):
        code = code[3:]
    
    # Remove file extensions
    if code.endswith('.txt') or code.endswith('.json'):
        code = code.rsplit('.', 1)[0]
    
    # Ensure 5 digits with leading zeros if it's numeric
    if code.isdigit():
        code = code.zfill(5)
    
    return code


def main():
    """Main function to load Excel file and download BRITE JSON files."""
    args = parse_args()
    
    # Set up paths
    script_dir = Path(__file__).parent
    excel_file = script_dir / "pombe_kegg_brite.xlsx"
    output_dir = args.output_dir if args.output_dir else script_dir / "brite_json"
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    
    # Determine which BRITE codes to download
    if args.codes:
        # User provided specific codes
        brite_codes = [normalize_brite_code(code) for code in args.codes]
        print(f"Downloading {len(brite_codes)} specific BRITE code(s): {', '.join(brite_codes)}\n")
    else:
        # Load all codes from Excel file
        if not excel_file.exists():
            print(f"Error: Excel file not found: {excel_file}")
            print("Please provide BRITE codes as arguments or ensure the Excel file exists.")
            sys.exit(1)
        
        brite_codes = get_brite_codes_from_excel(excel_file)
    
    # Download each BRITE JSON file
    successful = 0
    failed = 0
    
    for brite_code in brite_codes:
        result = download_brite_json(brite_code, output_dir, organism=args.organism)
        
        if result:
            successful += 1
        else:
            failed += 1
        
        # Be nice to the server
        time.sleep(1)
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Download complete!")
    print(f"  Successful: {successful}")
    print(f"  Failed: {failed}")
    print(f"  Output directory: {output_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
