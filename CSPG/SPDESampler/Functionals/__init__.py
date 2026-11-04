"""Quantities of interest for PDE solutions."""

from .abstract_qoi import AbstractQOI
from .average import AverageQOI
from .box_avg import BoxAverage

__all__ = ["AbstractQOI", "AverageQOI", "BoxAverage"]

