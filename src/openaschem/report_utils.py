"""Shared numerical serialization for offline reports."""

import math


def finite_number(value):
    """Return a finite float or None for JSON-safe display metadata."""
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None
