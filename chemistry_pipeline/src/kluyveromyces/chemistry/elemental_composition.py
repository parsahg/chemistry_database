from .formula_parser import parse_formula


# Standard atomic weights (g/mol)
ATOMIC_WEIGHTS = {
    "H": 1.008,   "He": 4.003,  "Li": 6.941,  "Be": 9.012,  "B": 10.811,
    "C": 12.011,  "N": 14.007,  "O": 15.999,  "F": 18.998,  "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.086, "P": 30.974,
    "S": 32.065,  "Cl": 35.453, "Ar": 39.948, "K": 39.098,  "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942,  "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.409,
    "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904,
    "Kr": 83.798, "Rb": 85.468, "Sr": 87.620, "Y": 88.906,  "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.950, "Tc": 98.000, "Ru": 101.070,"Rh": 102.910,
    "Pd": 106.420,"Ag": 107.870,"Cd": 112.410,"In": 114.820,"Sn": 118.710,
    "Sb": 121.760,"Te": 127.600,"I": 126.900, "Xe": 131.290,"Cs": 132.910,
    "Ba": 137.330,"La": 138.910,"Ce": 140.120,"Pr": 140.910,"Nd": 144.240,
    "Pm": 145.000,"Sm": 150.360,"Eu": 151.960,"Gd": 157.250,"Tb": 158.930,
    "Dy": 162.500,"Ho": 164.930,"Er": 167.260,"Tm": 168.930,"Yb": 173.050,
    "Lu": 174.970,"Hf": 178.490,"Ta": 180.950,"W": 183.840, "Re": 186.210,
    "Os": 190.230,"Ir": 192.220,"Pt": 195.080,"Au": 196.970,"Hg": 200.590,
    "Tl": 204.380,"Pb": 207.200,"Bi": 208.980,"Po": 209.000,"At": 210.000,
    "Rn": 222.000,"Fr": 223.000,"Ra": 226.000,"Ac": 227.000,"Th": 232.040,
    "Pa": 231.040,"U": 238.030, "Np": 237.000,"Pu": 244.000,"Am": 243.000,
    "Cm": 247.000,"Bk": 247.000,"Cf": 251.000,"Es": 252.000,"Fm": 257.000,
    "Md": 258.000,"No": 259.000,"Lr": 266.000,"Rf": 267.000,"Db": 268.000,
    "Sg": 269.000,"Bh": 270.000,"Hs": 277.000,"Mt": 278.000,"Ds": 281.000,
    "Rg": 282.000,"Cn": 285.000,"Nh": 286.000,"Fl": 289.000,"Mc": 290.000,
    "Lv": 293.000,"Ts": 294.000,"Og": 294.000,
}


def molecular_weight(composition: dict[str, int]) -> float:
    """
    Calculate molecular weight from an elemental composition.

    Parameters
    ----------
    composition : dict[str, int]
        Element counts, for example:
        {"C": 12, "H": 22, "O": 11}

    Returns
    -------
    float
        Molecular weight in g/mol.
    """

    molecular_mass = 0.0

    for element, count in composition.items():

        if element not in ATOMIC_WEIGHTS:
            raise ValueError(
                f"Atomic weight not available for element: {element}"
            )

        molecular_mass += ATOMIC_WEIGHTS[element] * count

    return molecular_mass


def molecular_weight_from_formula(formula: str) -> float:
    """
    Calculate molecular weight directly from a chemical formula.

    Example
    -------
    C12H22O11 -> 342.297 g/mol
    """

    composition = parse_formula(formula)

    return molecular_weight(composition)

def elemental_mass_fractions(
    composition: dict[str, int],
) -> dict[str, float]:
    """
    Calculate the mass fraction of each element in a compound.

    Parameters
    ----------
    composition : dict[str, int]
        Element counts, for example:
        {"C": 12, "H": 22, "O": 11}

    Returns
    -------
    dict[str, float]
        Elemental mass fractions.

    Example
    -------
    C12H22O11 ->
    {
        "C": ...,
        "H": ...,
        "O": ...
    }
    """

    total_mass = molecular_weight(composition)

    fractions = {}

    for element, count in composition.items():
        if element not in ATOMIC_WEIGHTS:
            raise ValueError(
                f"Atomic weight not available for element: {element}"
            )

        element_mass = count * ATOMIC_WEIGHTS[element]

        fractions[element] = element_mass / total_mass

    return fractions


def elemental_mass_fractions_from_formula(
    formula: str,
) -> dict[str, float]:
    """
    Calculate elemental mass fractions directly from a chemical formula.
    """

    composition = parse_formula(formula)

    return elemental_mass_fractions(composition)

def elemental_masses(
    composition: dict[str, int],
    material_mass: float,
) -> dict[str, float]:
    """
    Calculate the mass of each element contained in a given
    mass of a compound.

    Parameters
    ----------
    composition : dict[str, int]
        Element counts in the compound.

    material_mass : float
        Mass of the compound in grams.

    Returns
    -------
    dict[str, float]
        Mass of each element in grams.
    """

    if material_mass < 0:
        raise ValueError("Material mass cannot be negative.")

    fractions = elemental_mass_fractions(composition)

    return {
        element: material_mass * fraction
        for element, fraction in fractions.items()
    }


def elemental_masses_from_formula(
    formula: str,
    material_mass: float,
) -> dict[str, float]:
    """
    Calculate elemental masses directly from a chemical formula
    and the mass of the compound.

    Example
    -------
    40 g of C12H22O11 -> elemental masses in grams.
    """

    composition = parse_formula(formula)

    return elemental_masses(
        composition,
        material_mass,
    )









































