"""Regulatory annotation: what binds a promoter, and how well that is actually known.

The scientific core of this subpackage is :mod:`sarna.regulatory.evidence`, which separates three
claims a tool might make about a factor and a promoter — the sequence could bind it, the region is
regulatory somewhere, the factor was seen bound in this cell line — and refuses to let the first
masquerade as the third.
"""
