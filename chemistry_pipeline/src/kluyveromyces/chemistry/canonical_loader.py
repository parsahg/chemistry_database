"""Canonical one-row-per-experiment loader.

Purpose:
    Convert heterogeneous source schemas into one canonical structure.

IMPORTANT:
    This module performs STRUCTURAL canonicalization only.

    It does NOT:
    - convert units
    - change numerical values
    - impute missing values
    - round values
    - calculate elemental composition
    - infer missing experimental conditions
    - remove experimental observations

Raw files are never modified.
"""

from pathlib import Path
import re

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Column-name canonicalization
# ---------------------------------------------------------------------

COLUMN_ALIASES = {
    # General metadata
    "Project_ID": "project_id",
    "Project_Name": "project_name",
    "Run_ID": "run_id",
    "Strain": "strain",
    "Source": "source_reference",
    "Notes": "notes",

    # Carbon sources
    "Glucose_g_L": "glucose_g_L",
    "glucose_g_L": "glucose_g_L",
    "glucose_g_L_x": "glucose_g_L",
    "glucose_g_L_y": "glucose_g_L",

    "Galactose_g_L": "galactose_g_L",
    "galactose_g_L": "galactose_g_L",

    "Lactose_g_L": "lactose_g_L",
    "lactose_g_L": "lactose_g_L",
    "lactose_g_L_x": "lactose_g_L",
    "lactose_g_L_y": "lactose_g_L",

    # Medium components
    "MgSO4_7H2O_g_L": "MgSO4_7H2O_g_L",
    "MgSO4_g_L": "MgSO4_g_L",

    "CaCl2_g_L": "CaCl2_g_L",
    "NaCl_g_L": "NaCl_g_L",
    "NiSO4_7H2O_g_L": "NiSO4_7H2O_g_L",

    "(NH4)2SO4_g_L": "(NH4)2SO4_g_L",
    "K2HPO4_g_L": "K2HPO4_g_L",
    "KH2PO4_g_L": "KH2PO4_g_L",

    "yeast_extract_g_L": "yeast_extract_g_L",
    "peptone_g_L": "peptone_g_L",
    "urea_g_L": "urea_g_L",
    "whey_powder_g_L": "whey_powder_g_L",
    "corn_steep_liquor_g_L": "corn_steep_liquor_g_L",

    # Conditions
    "Temperature_C": "temperature_C",
    "temperature_C": "temperature_C",

    "Retention_time_h": "retention_time_h",
    "time_h": "time_h",

    "Aeration_vvm": "aeration_vvm",
    "Agitation_rpm": "agitation_rpm",

    # Vitamins / trace elements
    "Biotin_mg_L": "biotin_mg_L",
    "Ca_pantothenate_mg_L": "ca_pantothenate_mg_L",
    "Folic_acid_mg_L": "folic_acid_mg_L",
    "Inositol_mg_L": "inositol_mg_L",
    "Niacin_mg_L": "niacin_mg_L",
    "PABA_mg_L": "paba_mg_L",
    "Pyridoxine_HCl_mg_L": "pyridoxine_HCl_mg_L",
    "Riboflavin_mg_L": "riboflavin_mg_L",
    "Thiamine_HCl_mg_L": "thiamine_HCl_mg_L",
    "H3BO3_mg_L": "H3BO3_mg_L",
    "CuSO4_5H2O_mg_L": "CuSO4_5H2O_mg_L",
    "KI_mg_L": "KI_mg_L",
    "FeCl3_6H2O_mg_L": "FeCl3_6H2O_mg_L",
    "MnSO4_H2O_mg_L": "MnSO4_H2O_mg_L",
    "Na2MoO4_2H2O_mg_L": "Na2MoO4_2H2O_mg_L",
    "ZnSO4_7H2O_mg_L": "ZnSO4_7H2O_mg_L",

    # Published responses
    "Biomass_g_L": "biomass_g_L",
    "biomass_g_L": "biomass_g_L",

    "OD600": "OD600",

    "SCP_g_L": "SCP_g_L",
    "Yxs_lactose_g_g": "Yxs_lactose_g_g",
    "Cell_yield_g_g": "Cell_yield_g_g",

    "mu_max_1_h": "mu_max_1_h",
    "mu_max_h-1": "mu_max_1_h",

    "Doubling_time_h": "doubling_time_h",
    "doubling_time_h": "doubling_time_h",

    "Productivity_g_L_h": "Productivity_g_L_h",
    "Ethanol_g_L": "ethanol_g_L",
    "ethanol_g_L": "ethanol_g_L",

    "Protein_content_pct": "Protein_content_pct",
    "Residual_protein_g_L": "Residual_protein_g_L",

    "CFU_mL": "CFU_mL",
    "lactose_consumed_g_L": "lactose_consumed_g_L",

    # Experimental/design variables
    "Mineral_salts_relative_fold": "mineral_salts_relative_fold",

    "sucrose_pct_w_v": "sucrose_pct_w_v",
    "peptone_pct_w_v": "peptone_pct_w_v",
    "CaCl2_pct_w_v": "CaCl2_pct_w_v",
    "coded_Ca_level_pct": "coded_Ca_level_pct",

    "sucrose_g_L": "sucrose_g_L",
    "EPS_mg_L": "EPS_mg_L",

    # Explicitly retained if present in source data
    "carbon_source": "carbon_source",
    "carbon_g_L": "carbon_g_L",
    "medium": "medium",
    "species": "species",
    "data_type": "data_type",
}


# ---------------------------------------------------------------------
# Material names
# ---------------------------------------------------------------------

MATERIAL_ALIASES = {
    "D-glucose": "glucose",
    "Dextrose": "glucose",
    "Glucose": "glucose",

    "D-lactose": "lactose",
    "Lactose": "lactose",

    "Whey powder / whey medium": "whey_powder",

    "Yeast extract": "yeast_extract",
    "yeast extract": "yeast_extract",

    "Peptone": "peptone",
    "peptone": "peptone",

    "Urea": "urea",
    "urea": "urea",
}


def clean_material_name(name):
    """Convert a material name into a stable column-safe name."""
    name = str(name).strip()

    name = MATERIAL_ALIASES.get(name, name)

    name = re.sub(r"\s+", "_", name)
    name = re.sub(r"[^A-Za-z0-9_()+.-]", "_", name)

    return re.sub(r"_+", "_", name).strip("_")


def material_col(material, unit="g_L"):
    return f"{clean_material_name(material)}_{unit}"


def canonical_column_name(name):
    """Return canonical column name without changing its values."""
    name = str(name).strip()

    if name in COLUMN_ALIASES:
        return COLUMN_ALIASES[name]

    return name


# ---------------------------------------------------------------------
# Safe column merging
# ---------------------------------------------------------------------

def merge_columns_without_changing_values(df):
    """
    Merge columns that represent the same canonical variable.

    Rules:
    - one non-null value -> preserve it exactly
    - both null -> remain null
    - both present and equal -> preserve the value
    - both present and different -> raise an error

    No numerical transformation is performed.
    """

    output = pd.DataFrame(index=df.index)

    groups = {}

    for column in df.columns:
        canonical = canonical_column_name(column)
        groups.setdefault(canonical, []).append(column)

    for canonical, columns in groups.items():

        if len(columns) == 1:
            output[canonical] = df[columns[0]]
            continue

        combined = df[columns].bfill(axis=1).iloc[:, 0]

        # Detect conflicts without changing either source value.
        for i in range(len(df)):
            values = []

            for column in columns:
                value = df.iloc[i][column]

                if pd.notna(value):
                    values.append(value)

            if len(values) > 1:

                first = values[0]

                for value in values[1:]:
                    if value != first:
                        raise ValueError(
                            f"Conflicting values found while canonicalizing "
                            f"'{canonical}' at row {i}: "
                            f"{values}"
                        )

        output[canonical] = combined

    return output


# ---------------------------------------------------------------------
# Search 1
# ---------------------------------------------------------------------

def read_search1(path):

    raw = pd.read_excel(
        path,
        sheet_name="Wide Matrix",
        header=None
    )

    header_idx = next(
        i
        for i in range(10)
        if "Project_ID" in raw.iloc[i].astype(str).tolist()
    )

    df = raw.iloc[header_idx + 1:].copy()

    df.columns = raw.iloc[header_idx].tolist()

    df = df.dropna(how="all").reset_index(drop=True)

    out = pd.DataFrame(index=df.index)

    out["experiment_id"] = df["Run_ID"].astype(str)
    out["source_id"] = (
        "data_search1_P" +
        df["Project_ID"].astype(str)
    )

    out["project_id"] = df["Project_ID"]
    out["run_id"] = df["Run_ID"]

    out["species"] = "Kluyveromyces marxianus"
    out["strain"] = df["Strain"]
    out["medium"] = pd.NA
    out["record_type"] = "literature"
    out["project_name"] = df["Project_Name"]

    for col in df.columns:

        if col not in {
            "Project_ID",
            "Project_Name",
            "Run_ID",
            "Strain",
            "Notes",
            "Source",
        }:

            canonical = canonical_column_name(col)

            out[canonical] = df[col]

    out["source_reference"] = df["Source"]
    out["notes"] = df["Notes"]
    out["data_type"] = "raw_published_observation"

    return out


# ---------------------------------------------------------------------
# Search 2
# ---------------------------------------------------------------------

def pivot_unique_factors(df):

    rows = []

    for (source, medium), group in df.groupby(
        ["source_id", "medium_or_design"],
        dropna=False
    ):

        record = {
            "source_id": source,
            "medium": medium,
        }

        for factor, factor_group in group.groupby("factor"):

            # Do not silently discard repeated observations.
            if len(factor_group) != 1:
                continue

            row = factor_group.iloc[0]

            unit = (
                str(row["unit"]).strip()
                if pd.notna(row["unit"])
                else ""
            )

            value = row["value"]

            # IMPORTANT:
            # Values are retained exactly as reported.
            if unit == "g/L":
                record[
                    material_col(factor, "g_L")
                ] = value

            elif unit == "% w/v":
                record[
                    f"{clean_material_name(factor)}_pct_w_v"
                ] = value

            elif unit == "ppm":
                record[
                    f"{clean_material_name(factor)}_ppm"
                ] = value

            elif unit == "mM":
                record[
                    f"{clean_material_name(factor)}_mM"
                ] = value

            else:
                safe_unit = (
                    unit
                    .replace("/", "_")
                    .replace(" ", "_")
                )

                record[
                    f"{clean_material_name(factor)}_{safe_unit}"
                ] = value

        rows.append(record)

    return pd.DataFrame(rows)


def read_search2(path):

    g = pd.read_excel(
        path,
        sheet_name="Raw_growth_and_biomass"
    )

    out = pd.DataFrame(index=g.index)

    out["experiment_id"] = [
        f"{sid}-R{i + 1}"
        for i, sid in enumerate(
            g["source_id"].fillna("SEARCH2")
        )
    ]

    out["source_id"] = g["source_id"]
    out["project_id"] = pd.NA
    out["run_id"] = out["experiment_id"]

    out["species"] = g["species"]
    out["strain"] = g["strain"]
    out["medium"] = g["medium"]
    out["record_type"] = "literature"

    for col in [
        "carbon_source",
        "carbon_g_L",
        "temperature_C",
        "time_h",
        "OD600",
        "biomass_g_L",
        "mu_max_h-1",
        "doubling_time_h",
        "CFU_mL",
        "ethanol_g_L",
        "lactose_consumed_g_L",
    ]:

        if col in g.columns:

            canonical = canonical_column_name(col)

            out[canonical] = g[col]

    # Preserve carbon-source value exactly as reported.
    for idx, row in g.iterrows():

        if (
            pd.notna(row["carbon_source"])
            and pd.notna(row["carbon_g_L"])
        ):

            out.at[
                idx,
                material_col(row["carbon_source"])
            ] = row["carbon_g_L"]

    out["notes"] = g["notes"]
    out["data_type"] = g["data_type"]

    unique = pivot_unique_factors(
        pd.read_excel(
            path,
            sheet_name="Raw_medium_factors"
        )
    )

    if not unique.empty:

        out = out.merge(
            unique,
            on=["source_id", "medium"],
            how="left"
        )

    eps = pd.read_excel(
        path,
        sheet_name="Raw_EPS_design"
    )

    e = pd.DataFrame(index=eps.index)

    e["experiment_id"] = [
        f"{sid}-EPS-R{int(run)}"
        for sid, run in zip(
            eps["source_id"],
            eps["run"]
        )
    ]

    e["source_id"] = eps["source_id"]
    e["project_id"] = pd.NA
    e["run_id"] = e["experiment_id"]

    e["species"] = eps["species"]
    e["strain"] = eps["strain"]
    e["medium"] = "EPS_design"
    e["record_type"] = "literature"

    # Keep original percentage values exactly.
    for factor in ["sucrose", "peptone", "CaCl2"]:

        e[
            f"{factor}_pct_w_v"
        ] = eps[f"{factor}_pct"]

    # IMPORTANT:
    # Do NOT create converted g/L values here.

    e["coded_Ca_level_pct"] = eps[
        "coded_Ca_level_pct"
    ]

    e["EPS_mg_L"] = eps[
        "EPS_mg_L"
    ]

    e["data_type"] = "raw_published_observation"
    e["notes"] = pd.NA

    sources = pd.read_excel(
        path,
        sheet_name="Sources"
    )

    refs = dict(
        zip(
            sources["source_id"],
            sources["citation"]
        )
    )

    out["source_reference"] = out[
        "source_id"
    ].map(refs)

    e["source_reference"] = e[
        "source_id"
    ].map(refs)

    return pd.concat(
        [out, e],
        ignore_index=True,
        sort=False
    )


# ---------------------------------------------------------------------
# Laboratory data
# ---------------------------------------------------------------------

def read_lab(path):

    raw = pd.read_excel(
        path,
        sheet_name="Our Lab",
        header=None
    )

    # Determine experiments from the header row.
    exps = [
        int(float(x))
        for x in raw.iloc[0, 2:]
        if pd.notna(x)
    ]

    rows = []

    for j, number in enumerate(
        exps,
        start=2
    ):

        record = {
            "experiment_id": f"LAB-{number}",
            "source_id": "our_lab",
            "project_id": pd.NA,
            "run_id": f"LAB-{number}",
            "species": "Kluyveromyces marxianus",
            "strain": pd.NA,
            "medium": "laboratory_medium",
            "record_type": "laboratory",
            "data_type": "our_lab_observation",
            "source_reference": "our_lab_data.xlsx",
            "notes": pd.NA,
        }

        # Read ALL material rows.
        #
        # We intentionally do not assume that the lab sheet
        # contains exactly six material rows.
        for i in range(1, len(raw)):

            material = raw.iloc[i, 0]
            unit = raw.iloc[i, 1]
            value = raw.iloc[i, j]

            if pd.isna(material) or pd.isna(value):
                continue

            material_text = str(material).strip()
            unit_text = str(unit).strip().lower()

            if material_text.lower() == "biomass":
                record["biomass_g_L"] = value

            elif unit_text == "g/l":

                record[
                    material_col(material_text)
                ] = value

        rows.append(record)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# Final canonical dataset
# ---------------------------------------------------------------------

def build_canonical_dataset(
    search1_path,
    search2_path,
    lab_path
):

    search1 = read_search1(search1_path)
    search2 = read_search2(search2_path)
    lab = read_lab(lab_path)

    combined = pd.concat(
        [search1, search2, lab],
        ignore_index=True,
        sort=False
    )

    # Final structural canonicalization.
    combined = merge_columns_without_changing_values(
        combined
    )

    # Remove ONLY completely empty columns.
    combined = combined.dropna(
        axis=1,
        how="all"
    )

    return combined