#!/usr/bin/env python3
"""
Extract OD data from all sheets of all Excel files starting at row with 'Time' in column B or C.

Handles:
- Flexible 'Time' detection
- Skips layout sheets
- Fixes dtype issues
"""



#Usage : python OD_READINGS_CONVERTER_TO_TSV_FROM_EXCEL.py




import pandas as pd
import argparse
from pathlib import Path
import re

def clean_sheet_name(sheet_name):
    base = sheet_name.split("-")[0].strip()
    cleaned = re.sub(r"\s+", "_", base.upper())
    return cleaned

def extract_from_b39(df_raw, output_path):
    try:
        # Safety check: cell B39 = row 38, col 1
        b39 = df_raw.iat[38, 1] if df_raw.shape[0] > 38 and df_raw.shape[1] > 1 else ""
        b39_str = str(b39).strip().lower() if pd.notna(b39) else ""
        if b39_str != "time":
            print(f"⚠️  Warning: B39 in {output_path.name} does not contain 'Time' (found: '{b39}')")

        # Extract from B39 onward
        df = df_raw.iloc[38:, 1:]
        df.columns = df.iloc[0]
        df = df.drop(df.index[0])
        df.columns = [str(c) if pd.notna(c) else f"COL_{i}" for i, c in enumerate(df.columns)]

        # Identify time column
        time_col = df.columns[0]

        # Stop at first empty or "00:00:00"
        stop_mask = df[time_col].astype(str).apply(lambda t: str(t).strip() == "" or str(t).strip() == "00:00:00")
        if stop_mask.any():
            df = df.loc[:stop_mask.idxmax() - 1]

        # Drop temperature and unnamed columns
        df = df.loc[:, ~df.columns.astype(str).str.contains("T°|Unnamed|EMPTY", case=False, na=False)]

        # Clean all non-time columns: commas → dots
        for col in df.columns[1:]:
            try:
                df[col] = df[col].astype(str).str.replace(",", ".", regex=False)
            except Exception:
                print(f"⚠️  Skipping malformed column: {col}")

        # Replace time with simple TIME = 0, 1, 2...
        df.insert(0, "TIME", range(len(df)))
        df = df.drop(columns=[time_col])

        # Save to TSV
        df.to_csv(output_path, sep="\t", index=False)
        print(f"✔️ Saved: {output_path.name}")

    except Exception as e:
        print(f"❌ Error processing {output_path.name}: {e}")

def process_excel_file(file_path):
    file = Path(file_path)
    if not file.exists() or file.suffix.lower() != ".xlsx":
        print(f"❌ Error: File not found or not an .xlsx file: {file}")
        return

    try:
        xls = pd.ExcelFile(file, engine="openpyxl")
        for sheet_name in xls.sheet_names:
            df_raw = xls.parse(sheet_name=sheet_name, header=None, dtype=str)

            base_sheet = clean_sheet_name(sheet_name)
            output_name = f"{base_sheet}.tsv"
            output_path = file.parent / output_name

            extract_from_b39(df_raw, output_path)

    except Exception as e:
        print(f"❌ Error reading {file.name}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract OD data from Excel file (starting at B39).")
    parser.add_argument("excel_file", help="Path to the Excel .xlsx file to process")
    args = parser.parse_args()
    process_excel_file(args.excel_file)
