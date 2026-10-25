"""Shared defaults and validation for the standard reference RT window."""

import math

DEFAULT_RT_LEFT_LIMIT = 50
DEFAULT_RT_RIGHT_LIMIT = 300


def validate_rt_limits(left, right):
    """Return finite, ordered RT limits in seconds before any output is created."""
    left, right = float(left), float(right)
    if not (math.isfinite(left) and math.isfinite(right) and 0 <= left < right):
        raise ValueError("RT limits must be finite and satisfy 0 <= left < right (seconds).")
    return left, right
