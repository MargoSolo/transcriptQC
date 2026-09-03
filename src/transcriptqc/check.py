"""The checks. One transcript in, one Result out; files are just many of those."""
from __future__ import annotations
import csv, re
from dataclasses import dataclass, asdict, field
from pathlib import Path
from .mane import MANE, load_mane, split_acc

STATUS_HELP = {
    "MANE_SELECT": "matches the current MANE Select transcript — aligned with the MANE reference set",
    "MANE_PLUS_CLINICAL": "a MANE Plus Clinical transcript — acceptable when this isoform is the clinically relevant one; say why",
    "OLD_VERSION": "right transcript, outdated version — c. coordinates may differ; update the version and re-check the HGVS",
    "OLD_VERSION_PLUS_CLINICAL": "a MANE Plus Clinical transcript at an outdated version — update the version",
    "NEWER_VERSION": "a version newer than the MANE snapshot — update the MANE snapshot before deciding",
    "MANE_SELECT_CHANGED": "this used to be MANE Select; the current MANE Select is a different accession — check coordinates carefully",
    "NON_MANE": "not a MANE transcript for this gene — review the transcript choice against MANE Select / MANE Plus Clinical; a non-MANE transcript can be legitimate with a documented rationale",
    "GENE_MISMATCH": "the transcript belongs to a different gene than stated — likely a copy-paste error",
    "GENE_NOT_IN_MANE": "protein-coding gene without a MANE transcript yet — no MANE reference to compare against",
    "UNKNOWN_TRANSCRIPT": "accession not in the MANE snapshot and no gene given — cannot tell which gene it belongs to (add a gene column)",
    "UNKNOWN_GENE": "gene symbol not found in MANE — check the symbol (HGNC) or spelling",
    "UNPARSEABLE": "could not parse the transcript accession",
}
PROBLEM = {"OLD_VERSION", "OLD_VERSION_PLUS_CLINICAL", "MANE_SELECT_CHANGED", "NON_MANE", "GENE_MISMATCH", "UNKNOWN_TRANSCRIPT", "UNKNOWN_GENE", "UNPARSEABLE", "NEWER_VERSION"}
HGVS = re.compile(r"^\s*(?P<acc>[A-Z]{2,4}_?\d+(?:\.\d+)?)\s*(?:\((?P<gene>[A-Za-z0-9\-]+)\))?\s*:\s*(?P<var>[cgnmpr]\..+?)\s*$")


@dataclass
class Result:
    input: str
    gene: str | None
    variant: str | None
    status: str
    mane_select: str | None = None
    mane_select_ensembl: str | None = None
    plus_clinical: list[str] = field(default_factory=list)
    gene_of_transcript: str | None = None
    action: str = ""
    note: str = ""
    problem: bool = False

    def to_dict(self): return asdict(self)


def _result(inp, gene, var, status, e=None, m: MANE | None = None, **kw) -> Result:
    r = Result(inp, gene, var, status, **kw)
    sel = m.select(gene) if (m and gene) else (e if e is not None and e.status == "MANE Select" else None)
    if sel is None and e is not None and m: sel = m.select(e.symbol)
    if sel: r.mane_select, r.mane_select_ensembl = sel.refseq, sel.ensembl
    g = gene or (e.symbol if e else None)
    if g and m: r.plus_clinical = [x.refseq for x in m.plus_clinical(g)]
    r.action = ACTION.get(status, "")
    if r.plus_clinical and status in ("MANE_SELECT", "OLD_VERSION", "NON_MANE"): r.note = (r.note + " " if r.note else "") + f"gene has MANE Plus Clinical transcript(s): {', '.join(r.plus_clinical)} — check which isoform carries the variant"
    r.problem = status in PROBLEM
    return r


ACTION = {"MANE_SELECT": "none", "MANE_PLUS_CLINICAL": "state why the Plus Clinical isoform is used", "OLD_VERSION": "update transcript version; re-validate c. position",
          "OLD_VERSION_PLUS_CLINICAL": "update transcript version", "NEWER_VERSION": "refresh the MANE snapshot", "MANE_SELECT_CHANGED": "move to the current MANE Select; re-map the variant",
          "NON_MANE": "review transcript choice against MANE Select / MANE Plus Clinical; document the rationale if retaining a non-MANE transcript", "GENE_MISMATCH": "fix gene or transcript", "GENE_NOT_IN_MANE": "no MANE transcript available; follow gene-specific or laboratory transcript-selection policy",
          "UNKNOWN_TRANSCRIPT": "add the gene symbol", "UNKNOWN_GENE": "check the gene symbol", "UNPARSEABLE": "fix the accession"}


def check_transcript(transcript: str, gene: str | None = None, mane: MANE | None = None, variant: str | None = None) -> Result:
    m = mane or load_mane(); t = transcript.strip(); g = gene.strip().upper() if gene else None
    base, ver = split_acc(t)
    if not re.fullmatch(r"(NM_|NR_|XM_|XR_|ENST)\d+", base): return _result(t, g, variant, "UNPARSEABLE", m=m)
    is_ens = base.startswith("ENST"); e = (m.by_ensembl_base if is_ens else m.by_refseq_base).get(base)
    if e is None:
        if base in m.changed:
            c = m.changed[base]; gg = g or (c.get("symbol") or "").upper() or None
            return _result(t, gg, variant, "MANE_SELECT_CHANGED", m=m, gene_of_transcript=c.get("symbol"),
                           note=f"MANE Select changed in MANE v{c.get('since')}: now {c.get('current')} ({c.get('current_ensembl')}), previously {c.get('old')}; CDS affected: {c.get('affects_cds')}")
        if g:
            if g in m.by_gene:
                if m.select(g) is None and m.plus_clinical(g): return _result(t, g, variant, "NON_MANE", m=m)
                return _result(t, g, variant, "NON_MANE", m=m)
            if g in m.not_in_mane: return _result(t, g, variant, "GENE_NOT_IN_MANE", m=m)
            return _result(t, g, variant, "UNKNOWN_GENE", m=m)
        return _result(t, g, variant, "UNKNOWN_TRANSCRIPT", m=m)
    if g and g != e.symbol.upper():
        # the accession is MANE for another gene
        return _result(t, g, variant, "GENE_MISMATCH", e, m, gene_of_transcript=e.symbol, note=f"{t} is the MANE {e.status.replace('MANE ', '')} transcript of {e.symbol}, not {g}")
    cur_base, cur_ver = split_acc(e.ensembl if is_ens else e.refseq)
    plus = e.status == "MANE Plus Clinical"
    if ver is None: return _result(t, g, variant, "MANE_PLUS_CLINICAL" if plus else "MANE_SELECT", e, m, gene_of_transcript=e.symbol, note="no version given; current is " + (e.ensembl if is_ens else e.refseq))
    if ver == cur_ver: return _result(t, g, variant, "MANE_PLUS_CLINICAL" if plus else "MANE_SELECT", e, m, gene_of_transcript=e.symbol)
    if ver < cur_ver: return _result(t, g, variant, "OLD_VERSION_PLUS_CLINICAL" if plus else "OLD_VERSION", e, m, gene_of_transcript=e.symbol, note=f"current version is {e.ensembl if is_ens else e.refseq}")
    return _result(t, g, variant, "NEWER_VERSION", e, m, gene_of_transcript=e.symbol, note=f"snapshot has {e.ensembl if is_ens else e.refseq}")


def check_hgvs(hgvs: str, mane: MANE | None = None) -> Result:
    """'NM_007294.3:c.5266dupC' or 'NM_007294.3(BRCA1):c.5266dupC'."""
    mm = HGVS.match(hgvs)
    if not mm: return Result(hgvs, None, None, "UNPARSEABLE", action=ACTION["UNPARSEABLE"], problem=True)
    return check_transcript(mm.group("acc"), mm.group("gene"), mane, mm.group("var"))


def mane_for_gene(gene: str, mane: MANE | None = None) -> dict:
    m = mane or load_mane(); g = gene.strip().upper(); sel = m.select(g); plus = m.plus_clinical(g)
    return {"gene": g, "found": bool(sel or plus), "mane_select": sel.refseq if sel else None, "mane_select_ensembl": sel.ensembl if sel else None,
            "mane_select_protein": sel.refseq_prot if sel else None, "plus_clinical": [(x.refseq, x.ensembl) for x in plus], "not_in_mane": g in m.not_in_mane, "release": m.release}


def check_file(path, mane: MANE | None = None) -> list[Result]:
    """CSV/TSV with columns gene, transcript, variant (any subset; or a single column of HGVS strings). One Result per row."""
    m = mane or load_mane(); p = Path(path); txt = p.read_text(encoding="utf-8")
    dialect = "excel-tab" if p.suffix.lower() in (".tsv", ".txt") or "\t" in txt.splitlines()[0] else "excel"
    rows = list(csv.DictReader(txt.splitlines(), dialect=dialect)); out = []
    if not rows: return out
    cols = {c.lower().strip(): c for c in rows[0]}
    tcol = next((cols[k] for k in ("transcript", "refseq", "accession", "nm", "transcript_id") if k in cols), None)
    gcol = next((cols[k] for k in ("gene", "symbol", "gene_symbol") if k in cols), None)
    vcol = next((cols[k] for k in ("variant", "hgvs", "hgvsc", "cdna", "c.") if k in cols), None)
    for r in rows:
        t = (r.get(tcol) or "").strip() if tcol else ""; g = (r.get(gcol) or "").strip() if gcol else None; v = (r.get(vcol) or "").strip() if vcol else None
        if not t and v and ":" in v: out.append(check_hgvs(v, m)); continue
        if not t and v and ":" not in v and not gcol: out.append(Result(v, g, v, "UNPARSEABLE", action=ACTION["UNPARSEABLE"], problem=True)); continue
        if ":" in t and not v: out.append(check_hgvs(t, m)); continue
        out.append(check_transcript(t, g or None, m, v or None))
    return out
