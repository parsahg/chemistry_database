from pathlib import Path
import json


PURE_COMPOUND = "pure_compound"
MIXTURE = "mixture"
BIOLOGICAL_MATERIAL = "biological_material"
UNKNOWN = "unknown"


PROJECT_ROOT = Path(__file__).resolve().parents[4]

MATERIALS_PATH = (
    PROJECT_ROOT
    / "chemistry_pipeline"
    / "data"
    / "reference"
    / "materials.json"
)


def normalize_material_name(name: str) -> str:
    """
    Normalize a material name for database lookup.
    """

    if not isinstance(name, str):
        raise TypeError("Material name must be a string.")

    name = name.strip()

    if not name:
        raise ValueError("Material name cannot be empty.")

    return " ".join(name.split()).lower()


def load_material_database() -> dict:
    """
    Load the material reference database from JSON.
    """

    if not MATERIALS_PATH.exists():
        raise FileNotFoundError(
            f"Material database not found: {MATERIALS_PATH}"
        )

    with MATERIALS_PATH.open(
        "r",
        encoding="utf-8",
    ) as file:
        database = json.load(file)

    if "classes" not in database:
        raise ValueError(
            "Material database must contain 'classes'."
        )

    if "materials" not in database:
        raise ValueError(
            "Material database must contain 'materials'."
        )

    return database


def get_material_record(name: str) -> dict:
    """
    Retrieve a material record from the JSON database.

    Materials that are not present in the database are
    returned as UNKNOWN.
    """

    normalized_name = normalize_material_name(name)

    database = load_material_database()

    for material in database["materials"]:

        names = material.get("original_names", [])

        normalized_names = {
            normalize_material_name(material_name)
            for material_name in names
        }

        canonical_name = material.get(
            "canonical_name",
            "",
        )

        if canonical_name:
            normalized_names.add(
                normalize_material_name(canonical_name)
            )

        if normalized_name in normalized_names:
            return {
                "name": name,
                "normalized_name": normalized_name,
                **material,
            }

    return {
        "name": name,
        "normalized_name": normalized_name,
        "class": UNKNOWN,
    }


def classify_material(name: str) -> str:
    """
    Return the classification of a material.
    """

    return get_material_record(name)["class"]


def create_material_record(
    name: str,
    material_class: str,
    source: str | None = None,
    source_reference: str | None = None,
) -> dict:
    """
    Create a standardized material record.
    """

    database = load_material_database()

    if material_class not in database["classes"]:
        raise ValueError(
            f"Unknown material class: {material_class}"
        )

    normalized_name = normalize_material_name(name)

    provenance = []

    if source is not None or source_reference is not None:
        provenance.append(
            {
                "source": source,
                "reference": source_reference,
            }
        )

    return {
        "material_id": None,
        "original_names": [name],
        "canonical_name": normalized_name,
        "class": material_class,
        "identifiers": {},
        "chemistry": {
            "formula": None,
            "molecular_weight": None,
        },
        "provenance": provenance,
        "confidence": "low",
        "notes": "",
    }

def update_material_record(
    material_id: str,
    updates: dict,
) -> dict:
    """
    Update selected fields of an existing material record.

    Existing fields are preserved unless explicitly updated.
    """

    database = load_material_database()

    for material in database["materials"]:
        if material.get("material_id") == material_id:

            material.update(updates)

            with MATERIALS_PATH.open(
                "w",
                encoding="utf-8",
            ) as file:
                json.dump(
                    database,
                    file,
                    ensure_ascii=False,
                    indent=2,
                )

            return material

    raise ValueError(
        f"Material ID not found: {material_id}"
    )

def generate_material_id(database: dict) -> str:
    """
    Generate the next available material ID.
    """

    existing_ids = []

    for material in database["materials"]:
        material_id = material.get("material_id")

        if (
            isinstance(material_id, str)
            and material_id.startswith("MAT-")
        ):
            try:
                number = int(material_id.split("-")[1])
                existing_ids.append(number)
            except ValueError:
                continue

    next_number = (
        max(existing_ids) + 1
        if existing_ids
        else 1
    )

    return f"MAT-{next_number:04d}"

def add_material_record(record: dict) -> dict:
    """
    Add a new material record to the database.
    """

    database = load_material_database()

    material_id = record.get("material_id")

    if not material_id:
        material_id = generate_material_id(database)
        record["material_id"] = material_id

    existing_ids = {
        material.get("material_id")
        for material in database["materials"]
    }

    if material_id in existing_ids:
        raise ValueError(
            f"Material ID already exists: {material_id}"
        )

    database["materials"].append(record)

    with MATERIALS_PATH.open(
        "w",
        encoding="utf-8",
    ) as file:
        json.dump(
            database,
            file,
            ensure_ascii=False,
            indent=2,
        )

    return record






























