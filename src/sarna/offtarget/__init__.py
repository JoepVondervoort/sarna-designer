"""Off-target screening: what else in the human genome or transcriptome a duplex could bind.

The engine designs against one promoter. Nothing in the design rules asks whether the resulting
duplex also matches somewhere else, and until this package existed the answer was simply absent —
every order export said ``off-target-screening: not-performed`` (ADR 0017).

Screening is **opt-in and asynchronous**: a run completes first, the researcher reads the designs,
and starts a search when they want one. It **never filters**. A finding is reported beside the
candidate it belongs to, with several separate measures rather than one score, because a single
perfect match elsewhere is a different experimental problem from forty weak partial ones and any
composite can rank them the same (ADR 0006, ADR 0019).
"""
