"""
CarrierCapture: Carrier capture and non-radiative recombination rate calculations.

A modern Python package for computing carrier capture rates in semiconductors
using multiphonon theory.
"""

from .__version__ import __version__
from .core import Potential, ConfigCoordinate, TransferCoordinate

__all__ = [
    "__version__",
    "Potential",
    "ConfigCoordinate",
    "TransferCoordinate",
]
