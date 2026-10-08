import re


def parse_formula(formula: str) -> dict[str, int]:
    """
    Parse a chemical molecular formula into element counts.

    Example
    -------
    C12H22O11 -> {"C": 12, "H": 22, "O": 11}
    """

    if not isinstance(formula, str):
        raise TypeError("Formula must be a string.")

    formula = formula.strip()

    if not formula:
        raise ValueError("Formula cannot be empty.")

    pattern = r"([A-Z][a-z]?)(\d*)"
    matches = re.findall(pattern, formula)

    if not matches:
        raise ValueError(f"Invalid chemical formula: {formula}")

    # Reconstruct the formula from the parsed tokens.
    # This prevents silently accepting malformed strings.
    reconstructed = "".join(
        element + count
        for element, count in matches
    )

    if reconstructed != formula:
        raise ValueError(f"Invalid chemical formula: {formula}")

    composition = {}

    for element, count in matches:
        count = int(count) if count else 1

        if element in composition:
            composition[element] += count
        else:
            composition[element] = count

    return composition