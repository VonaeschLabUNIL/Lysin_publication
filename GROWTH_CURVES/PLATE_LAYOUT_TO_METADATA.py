# PLATE_LAYOUT_TO_METADATA.py

import pandas as pd
import os
import sys
import argparse


# python PLATE_LAYOUT_TO_METADATA.py Plate_Layouts_Rotated.xlsx --output metadata.tsv



def extract_plate_metadata(sheet_df, plate_number):
    metadata = []

    for row in range(1, 13):  # Rows 1-12 (actual wells)
        well_number = sheet_df.iloc[row, 0]
        for col in range(1, 9):  # Columns B-I (1-8, letters A-H)
            well_letter = sheet_df.iloc[0, col]
            content = sheet_df.iloc[row, col]
            well = f"{well_letter}{int(well_number)}"

            if isinstance(content, str) and content.strip().upper() not in ["EMPTY", "BLANK"]:
                parts = content.split("-")
                labid     = parts[0].strip() if len(parts) > 0 else "EMPTY"
                species   = parts[1].strip() if len(parts) > 1 else "EMPTY"
                culture   = parts[2].strip() if len(parts) > 2 else "EMPTY"
                treatment = parts[3].strip() if len(parts) > 3 else "EMPTY"
                medium    = parts[4].strip() if len(parts) > 4 else "EMPTY"
                rep       = parts[5].strip() if len(parts) > 5 else "EMPTY"
                empty     = "USED"
            else:
                labid = species = culture = treatment = medium = rep = content.strip().upper() if isinstance(content, str) else "EMPTY"
                empty = labid  # Either "EMPTY" or "BLANK"

            metadata.append({
                "PLATE": plate_number,
                "PLATE_NAME": f"PLATE_{plate_number}",
                "WELL": well,
                "WELL_LETTER": well_letter,
                "WELL_NUMBER": int(well_number),
                "EMPTY": empty,
                "LABID": labid,
                "SPECIES": species,
                "CULTURE": culture,
                "TREATMENT": treatment,
                "MEDIUM": medium,
                "REP": rep
            })
    return metadata

def convert_excel_to_metadata(excel_file, output_tsv):
    xls = pd.ExcelFile(excel_file)
    all_metadata = []

    for i, sheet_name in enumerate(xls.sheet_names, start=1):
        df = xls.parse(sheet_name, header=None)
        all_metadata.extend(extract_plate_metadata(df, i))

    metadata_df = pd.DataFrame(all_metadata)
    metadata_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Metadata saved to: {output_tsv}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert rotated 96-well plate layouts (dash-separated metadata) to TSV.",
        epilog="Example: python PLATE_LAYOUT_TO_METADATA.py Plate_Layouts_Rotated.xlsx --output metadata.tsv"
    )
    parser.add_argument("excel_file", help="Path to the Excel file containing plate layouts")
    parser.add_argument("--output", default="metadata_output.tsv", help="Output TSV file path (default: metadata_output.tsv)")
    args = parser.parse_args()

    convert_excel_to_metadata(args.excel_file, args.output)

if __name__ == "__main__":
    main()
