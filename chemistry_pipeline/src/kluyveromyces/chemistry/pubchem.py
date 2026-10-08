from pathlib import Path

import pandas as pd
import requests

from .elemental_composition import elemental_masses
from .formula_parser import parse_formula

from .material_database import (
    PURE_COMPOUND,
    MIXTURE,
    BIOLOGICAL_MATERIAL,
    UNKNOWN,
    get_material_record,
    normalize_material_name,
)

PUBCHEM_BASE_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
)

PROJECT_ROOT = Path(__file__).resolve().parents[4]

PROJECT_ROOT = Path(__file__).resolve().parents[4]

CACHE_PATH = (
    PROJECT_ROOT
    / "chemistry_pipeline"
    / "data"
    / "reference"
    / "pubchem_cache.csv"
)

def normalize_compound_name(name: str) -> str:
    """
    Normalize a compound name for lookup and caching.

    This only removes harmless formatting differences.
    It does not perform chemical synonym resolution.
    """

    if not isinstance(name, str):
        raise TypeError("Compound name must be a string.")

    name = name.strip()

    if not name:
        raise ValueError("Compound name cannot be empty.")

    # Normalize repeated whitespace
    name = " ".join(name.split())

    # Normalize case for lookup
    return name.lower()

def get_compound(name: str) -> dict:
    """
    Retrieve basic compound information from PubChem by name.
    """

    if not isinstance(name, str):
        raise TypeError("Compound name must be a string.")

    name = name.strip()

    if not name:
        raise ValueError("Compound name cannot be empty.")

    url = (
        f"{PUBCHEM_BASE_URL}/compound/name/"
        f"{requests.utils.quote(name)}/property/"
        "MolecularFormula,MolecularWeight/JSON"
    )

    response = requests.get(url, timeout=15)

    if response.status_code == 404:
        raise ValueError(
            f"Compound not found in PubChem: {name}"
        )

    response.raise_for_status()

    data = response.json()

    properties = data["PropertyTable"]["Properties"][0]

    return {
        "name": name,
        "cid": properties["CID"],
        "formula": properties["MolecularFormula"],
        "molecular_weight": properties["MolecularWeight"],
    }

def get_synonyms(cid: int) -> list[str]:
    """
    Retrieve synonyms for a PubChem compound by CID.
    """

    url = (
        f"{PUBCHEM_BASE_URL}/compound/cid/"
        f"{cid}/synonyms/JSON"
    )

    response = requests.get(url, timeout=15)

    response.raise_for_status()

    data = response.json()

    return data["InformationList"]["Information"][0]["Synonym"]

def resolve_compound_name(name: str) -> dict:
    """
    Resolve a compound name through PubChem.

    Returns the canonical compound information together
    with the original input name.
    """

    normalized_name = normalize_compound_name(name)

    compound = get_compound(normalized_name)

    synonyms = get_synonyms(compound["cid"])

    normalized_synonyms = {
        normalize_compound_name(synonym)
        for synonym in synonyms
    }

    if normalized_name not in normalized_synonyms:
        raise ValueError(
            f"PubChem did not confirm '{name}' "
            f"as a synonym of CID {compound['cid']}."
        )

    return {
        "input_name": name,
        "canonical_name": compound["name"],
        "cid": compound["cid"],
        "formula": compound["formula"],
        "molecular_weight": compound["molecular_weight"],
    }

def load_cache() -> pd.DataFrame:
    """
    Load the local PubChem cache.

    Creates an empty cache if the file does not exist
    or is empty.
    """

    columns = [
        "name",
        "cid",
        "formula",
        "molecular_weight",
    ]

    if not CACHE_PATH.exists():
        return pd.DataFrame(columns=columns)

    if CACHE_PATH.stat().st_size == 0:
        return pd.DataFrame(columns=columns)

    return pd.read_csv(CACHE_PATH)


def save_to_cache(compound: dict) -> None:
    """
    Add a compound to the local PubChem cache.
    """

    cache = load_cache()

    # Remove an existing record for the same name.
    cache = cache[
        cache["name"].str.lower() != compound["name"].lower()
    ]

    cache = pd.concat(
        [cache, pd.DataFrame([compound])],
        ignore_index=True,
    )

    CACHE_PATH.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    cache.to_csv(
        CACHE_PATH,
        index=False,
    )

def get_compound_cached(name: str) -> dict:
    """
    Retrieve compound information using the local cache.

    If the compound is already cached, no PubChem request
    is made. Otherwise, the compound is resolved through
    PubChem and the result is stored locally.
    """

    normalized_name = normalize_compound_name(name)

    cache = load_cache()

    matches = cache[
        cache["name"].apply(normalize_compound_name)
        == normalized_name
    ]

    if not matches.empty:
        row = matches.iloc[0]

        return {
            "name": row["name"],
            "cid": int(row["cid"]),
            "formula": row["formula"],
            "molecular_weight": float(row["molecular_weight"]),
        }

    resolved = resolve_compound_name(name)

    compound = {
        "name": resolved["canonical_name"],
        "cid": resolved["cid"],
        "formula": resolved["formula"],
        "molecular_weight": resolved["molecular_weight"],
    }

    save_to_cache(compound)

    return compound

def compound_elemental_masses(
    name: str,
    material_mass: float,
) -> dict:
    """
    Retrieve a compound and calculate its elemental masses.
    """

    compound = get_compound_cached(name)

    composition = parse_formula(
        compound["formula"]
    )

    elemental_masses_result = elemental_masses(
        composition,
        material_mass,
    )

    return {
        "compound": compound,
        "elemental_masses": elemental_masses_result,
    }

def resolve_material(name: str) -> dict:
    """
    Resolve a material into a classification and, when possible,
    chemical information.

    Known mixtures and biological materials are handled locally.
    Unknown materials are checked against PubChem.
    """

    material = get_material_record(name)

    material_class = material["class"]

    if material_class in {
        MIXTURE,
        BIOLOGICAL_MATERIAL,
    }:
        return {
            "name": name,
            "normalized_name": material["normalized_name"],
            "class": material_class,
            "compound": None,
        }

    try:
        compound = get_compound_cached(name)

    except (ValueError, requests.RequestException):
        return {
            "name": name,
            "normalized_name": material["normalized_name"],
            "class": UNKNOWN,
            "compound": None,
        }

    return {
        "name": name,
        "normalized_name": material["normalized_name"],
        "class": PURE_COMPOUND,
        "compound": compound,
    }

def enrich_material_from_pubchem(material_id: str) -> dict:
    """
    Enrich a pure-compound material record using PubChem.

    The material must already exist in materials.json and
    must be classified as a pure compound.
    """

    from .material_database import (
        load_material_database,
        update_material_record,
    )

    database = load_material_database()

    material = None

    for record in database["materials"]:
        if record.get("material_id") == material_id:
            material = record
            break

    if material is None:
        raise ValueError(
            f"Material ID not found: {material_id}"
        )

    if material["class"] != PURE_COMPOUND:
        raise ValueError(
            f"Material '{material_id}' is not a pure compound."
        )

    name = material["canonical_name"]

    compound = get_compound_cached(name)

    updates = {
        "identifiers": {
            **material.get("identifiers", {}),
            "pubchem_cid": compound["cid"],
        },
        "chemistry": {
            **material.get("chemistry", {}),
            "formula": compound["formula"],
            "molecular_weight": compound["molecular_weight"],
        },
    }

    return update_material_record(
        material_id,
        updates,
    )

def enrich_all_pure_compounds() -> list[dict]:
    """
    Enrich all pure compounds in the material database
    using PubChem.

    Mixtures, biological materials, and unknown materials
    are not sent to PubChem.
    """

    from .material_database import load_material_database

    database = load_material_database()

    results = []

    for material in database["materials"]:

        if material["class"] != PURE_COMPOUND:
            continue

        material_id = material["material_id"]

        try:
            enriched = enrich_material_from_pubchem(
                material_id
            )

            results.append({
                "material_id": material_id,
                "name": enriched["canonical_name"],
                "status": "enriched",
                "cid": enriched["identifiers"].get(
                    "pubchem_cid"
                ),
                "formula": enriched["chemistry"].get(
                    "formula"
                ),
                "molecular_weight": enriched["chemistry"].get(
                    "molecular_weight"
                ),
            })

        except (ValueError, requests.RequestException) as error:

            results.append({
                "material_id": material_id,
                "name": material["canonical_name"],
                "status": "failed",
                "error": str(error),
            })

    return results

def resolve_and_enrich_material(name: str) -> dict:
    """
    Resolve a material using the local material database
    and PubChem when appropriate.

    Known mixtures and biological materials remain protected
    from molecular-formula lookup.
    """

    material = get_material_record(name)

    material_class = material["class"]

    # Known non-pure materials must not be treated
    # as single molecular compounds.
    if material_class in {
        MIXTURE,
        BIOLOGICAL_MATERIAL,
    }:
        return {
            "name": name,
            "class": material_class,
            "status": "requires_source_composition",
            "compound": None,
        }

    # Already known pure compound
    if material_class == PURE_COMPOUND:
        compound = get_compound_cached(
            material["canonical_name"]
        )

        return {
            "name": name,
            "class": PURE_COMPOUND,
            "status": "resolved",
            "compound": compound,
        }

    # Unknown material: try PubChem
    try:
        compound = get_compound_cached(name)

    except (ValueError, requests.RequestException):
        return {
            "name": name,
            "class": UNKNOWN,
            "status": "unresolved",
            "compound": None,
        }

    return {
        "name": name,
        "class": PURE_COMPOUND,
        "status": "resolved",
        "compound": compound,
    }

def resolve_and_register_material(name: str) -> dict:
    """
    Resolve a material and register it in the local
    material database when it is identified as a
    pure compound.
    """

    from .material_database import (
        add_material_record,
        load_material_database,
    )

    # First check the local database.
    material = get_material_record(name)

    if material["class"] != UNKNOWN:
        return resolve_and_enrich_material(name)

    # Try PubChem.
    try:
        compound = get_compound_cached(name)

    except (ValueError, requests.RequestException):
        return {
            "name": name,
            "class": UNKNOWN,
            "status": "unresolved",
            "compound": None,
        }

    # Check whether another material already contains
    # this compound as a synonym/canonical name.
    database = load_material_database()

    normalized_name = normalize_material_name(name)

    for existing in database["materials"]:

        names = existing.get("original_names", [])

        normalized_names = {
            normalize_material_name(value)
            for value in names
        }

        canonical_name = existing.get(
            "canonical_name"
        )

        if canonical_name:
            normalized_names.add(
                normalize_material_name(canonical_name)
            )

        if normalized_name in normalized_names:
            return resolve_and_enrich_material(
                name
            )

    record = {
        "material_id": None,
        "original_names": [name],
        "canonical_name": compound["name"],
        "class": PURE_COMPOUND,
        "identifiers": {
            "pubchem_cid": compound["cid"],
        },
        "chemistry": {
            "formula": compound["formula"],
            "molecular_weight": compound[
                "molecular_weight"
            ],
        },
        "provenance": [
            {
                "source": "PubChem",
                "reference": (
                    f"CID {compound['cid']}"
                ),
            }
        ],
        "confidence": "high",
        "notes": "",
    }

    record = add_material_record(record)

    return {
        "name": name,
        "class": PURE_COMPOUND,
        "status": "registered",
        "compound": compound,
        "material": record,
    }
































































