# Changelog

## 0.1.2 (2026-09-03)

'action' renamed to 'suggested check' everywhere (text, CSV, HTML; JSON key `suggested_check`): the tool flags what needs review, it does not prescribe clinical action.

## 0.1.1 (2026-09-03)

Wording: statuses describe alignment with the MANE reference set rather than 'right for clinical reporting'; NON_MANE asks for review and a documented rationale; GENE_NOT_IN_MANE defers to gene-specific or laboratory policy; accurate GRCh37 statement. PyPI release; Zenodo archiving enabled.

## 0.1.0 (2026-09-03)

First release: MANE v1.5 bundled; transcript / HGVS / gene / file checks; 12 statuses; text, CSV, JSON, HTML reports; --fail gate.
