
"""
Utility helpers shared by plotting scripts.

Currently provides:
- parse_index_list:  "1,3,5-10" -> [1,3,5,6,7,8,9,10]
- parse_atom_indices: ["1-3,5", "7-8"] -> [1,2,3,5,7,8]
- compute_point_density: per-point 2D density for density-colored scatter
"""

# Above this many points, skip the O(N^2) Gaussian KDE and use the
# histogram fallback instead, to keep density coloring responsive.
_KDE_MAX_POINTS = 10_000


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


def compute_point_density(x, y, *, bins: int = 64):
    """
    Estimate a per-point 2D density for density-colored scatter plots.

    Returns a float array of the same length as `x`/`y`; larger values mark
    points sitting in denser regions. A Gaussian KDE is used when SciPy is
    available and the point count is modest; otherwise (or for very large
    inputs, or a singular KDE covariance) a 2D-histogram bin-count fallback
    is used, so SciPy is not a hard dependency.

    # Arguments
    - `x`, `y`: equal-length sequences of point coordinates.
    - `bins`: number of bins per axis for the histogram fallback.

    # Returns
    - `numpy.ndarray` of densities, aligned with the input order.
    """
    import numpy as np

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = x.size
    if n == 0:
        return np.zeros(0, dtype=float)
    if n == 1:
        return np.ones(1, dtype=float)

    if n <= _KDE_MAX_POINTS:
        try:
            from scipy.stats import gaussian_kde

            xy = np.vstack([x, y])
            return np.asarray(gaussian_kde(xy)(xy), dtype=float)
        except (ImportError, np.linalg.LinAlgError, ValueError):
            # SciPy missing, or singular covariance (e.g. all points
            # collinear / identical): fall through to the histogram path.
            pass

    counts, xedges, yedges = np.histogram2d(x, y, bins=bins)
    ix = np.clip(np.searchsorted(xedges, x, side="right") - 1, 0, counts.shape[0] - 1)
    iy = np.clip(np.searchsorted(yedges, y, side="right") - 1, 0, counts.shape[1] - 1)
    return counts[ix, iy].astype(float)

