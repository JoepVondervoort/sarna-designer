"""Shared plumbing for talking to external data sources.

Four modules in this package reach the network (Ensembl, NCBI BLAST, the install check and
now GTEx) and each needs the same three behaviours: cache every response on disk so a result
is reproducible without the network, retry what is worth retrying, and refuse outright when
offline rather than silently returning nothing. :mod:`sarna.sources.json_http` is that
behaviour, written once.

:class:`sarna.genome.ensembl.EnsemblClient` predates this module and still carries its own
copy. Migrating it is a separate change: its retry semantics are load-bearing (Ensembl's
public endpoint serves an HTML error page with a 500 under load, and the comment explaining
that is the reason the retry loop has the shape it has) so it should move only with its own
test run, not as a side effect of adding a new source.
"""
