
"""
Utility helpers shared by plotting scripts.

Currently provides:
- parse_index_list:  "1,3,5-10" -> [1,3,5,6,7,8,9,10]
- parse_atom_indices: ["1-3,5", "7-8"] -> [1,2,3,5,7,8]
"""


def parse_index_list(spec: str) -> list[int]:
    """
    Parse a string like "1,3,5-10" into a list of integer indices.
    Supports comma-separated integers and ranges with hyphens.
    """
    indices: list[int] = []
    for token in spec.split(","):
        t = token.strip()
        if not t:
            continue
        if "-" in t:
            parts = token.split("-", 1)
            if len(parts) != 2:
                raise ValueError(f"Invalid range token: {t}")
            try:
                start_idx = int(parts[0].strip())
                end_idx = int(parts[1].strip())
            except ValueError as e:
                raise ValueError(f"Invalid integer in range: {t}") from e
            if start_idx > end_idx:
                raise ValueError(f"Range start greater than end: {t}")
            indices.extend(range(start_idx, end_idx + 1))
        else:
            try:
                indices.append(int(t))
            except ValueError as e:
                raise ValueError(f"Invalid integer token: {t}") from e
    return indices


def parse_atom_indices(atoms_args: list[str]) -> list[int]:
    """
    Parse atom-indices argument that may include ranges like '1-5'
    and comma-separated tokens, e.g. ['1-3,5,7-8'].
    """
    indices: list[int] = []
    for raw_token in atoms_args:
        if not raw_token:
            continue
        for token in raw_token.split(","):
            t = token.strip()
            if not t:
                continue
            if "-" in t:
                parts = t.split("-", 1)
                if len(parts) != 2:
                    raise ValueError(f"Invalid range token: {t}")
                try:
                    start_idx = int(parts[0].strip())
                    end_idx = int(parts[1].strip())
                except ValueError as e:
                    raise ValueError(f"Invalid integer in range: {t}") from e
                if start_idx > end_idx:
                    raise ValueError(f"Range start greater than end: {t}")
                indices.extend(range(start_idx, end_idx + 1))
            else:
                try:
                    indices.append(int(t))
                except ValueError as e:
                    raise ValueError(f"Invalid integer token: {t}") from e
    return indices

