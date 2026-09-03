"""The MANE snapshot: bundled (versioned inside the package) or a newer release downloaded to a cache."""
from __future__ import annotations
import gzip, os, re, urllib.request
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path

FTP = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human"
BUNDLED = "1.5"


def split_acc(acc: str) -> tuple[str, int | None]:
    """'NM_007294.3' → ('NM_007294', 3); 'ENST00000357654' → ('ENST00000357654', None)."""
    acc = acc.strip()
    m = re.fullmatch(r"([A-Z]{2,4}_?\d+)(?:\.(\d+))?", acc)
    if not m: return acc, None
    return m.group(1), int(m.group(2)) if m.group(2) else None


@dataclass
class Entry:
    symbol: str
    gene_id: str
    hgnc: str
    name: str
    refseq: str
    refseq_prot: str
    ensembl: str
    ensembl_prot: str
    status: str            # "MANE Select" | "MANE Plus Clinical"
    chrom: str
    start: int
    end: int
    strand: str


@dataclass
class MANE:
    release: str
    entries: list[Entry]
    changed: dict[str, dict] = field(default_factory=dict)         # old RefSeq/Ensembl base accession → {current, old, since, affects_cds} (from changed_select_accessions)
    not_in_mane: set[str] = field(default_factory=set)             # protein-coding gene symbols without a MANE transcript
    by_gene: dict[str, list[Entry]] = field(default_factory=dict)
    by_refseq_base: dict[str, Entry] = field(default_factory=dict)
    by_ensembl_base: dict[str, Entry] = field(default_factory=dict)

    def __post_init__(self):
        for e in self.entries:
            self.by_gene.setdefault(e.symbol.upper(), []).append(e)
            self.by_refseq_base[split_acc(e.refseq)[0]] = e; self.by_ensembl_base[split_acc(e.ensembl)[0]] = e

    def select(self, gene: str) -> Entry | None:
        return next((e for e in self.by_gene.get(gene.upper(), []) if e.status == "MANE Select"), None)

    def plus_clinical(self, gene: str) -> list[Entry]:
        return [e for e in self.by_gene.get(gene.upper(), []) if e.status == "MANE Plus Clinical"]


def _open(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def parse_summary(path: Path, release: str) -> MANE:
    entries = []
    with _open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"): header = line.lstrip("#").rstrip("\n").split("\t"); continue
            f = dict(zip(header, line.rstrip("\n").split("\t")))
            entries.append(Entry(f["symbol"], f["NCBI_GeneID"], f["HGNC_ID"], f["name"], f["RefSeq_nuc"], f["RefSeq_prot"], f["Ensembl_nuc"], f["Ensembl_prot"], f["MANE_status"],
                                 f["GRCh38_chr"], int(f["chr_start"]), int(f["chr_end"]), f["chr_strand"]))
    return MANE(release, entries)


def _add_extras(m: MANE, changed_path: Path | None, notin_path: Path | None):
    if changed_path and changed_path.exists():
        with _open(changed_path) as fh:
            hdr = None
            for line in fh:
                if line.startswith("#"): hdr = line.lstrip("#").rstrip("\n").split("\t"); continue
                f = line.rstrip("\n").split("\t")
                if hdr and len(f) >= len(hdr):
                    d = dict(zip(hdr, f)); rec = {"symbol": d.get("Symbol"), "current": d.get("Current_MANE_Select_RefSeq"), "current_ensembl": d.get("Current_MANE_Select_Ensembl"),
                                                  "old": d.get("Old_MANE_Select_RefSeq"), "old_ensembl": d.get("Old_MANE_Select_Ensembl"), "since": d.get("Current_MANE_Version"), "affects_cds": d.get("Update_Affects_CDS")}
                    for k in ("old", "old_ensembl"):
                        if rec.get(k): m.changed[split_acc(rec[k])[0]] = rec
    if notin_path and notin_path.exists():
        with _open(notin_path) as fh:
            for line in fh:
                if line.startswith("#"): continue
                f = line.rstrip("\n").split("\t")   # GeneID · HGNC_id · gene_symbol · status
                if len(f) >= 3 and f[2]: m.not_in_mane.add(f[2].upper())


def load_mane(release: str = "bundled") -> MANE:
    """`bundled` = the snapshot shipped with the package (v1.5). Any other tag ('1.5', '1.6', 'current') downloads that release's
    summary into ~/.cache/transcriptqc once; results then depend on that file, which is recorded in every report."""
    if release == "bundled":
        d = resources.files("transcriptqc") / "data"
        m = parse_summary(Path(str(d / f"MANE.GRCh38.v{BUNDLED}.summary.txt.gz")), f"v{BUNDLED} (bundled)")
        _add_extras(m, Path(str(d / f"MANE.GRCh38.v{BUNDLED}.changed_select_accessions.txt.gz")), Path(str(d / f"MANE.GRCh38.v{BUNDLED}.protein_coding_genes_not_in_mane.txt.gz"))); return m
    cache = Path(os.environ.get("TRANSCRIPTQC_CACHE", Path.home() / ".cache" / "transcriptqc")); cache.mkdir(parents=True, exist_ok=True)
    sub = "current" if release == "current" else f"release_{release}"
    listing = urllib.request.urlopen(f"{FTP}/{sub}/", timeout=60).read().decode()
    fn = re.search(r"MANE\.GRCh38\.v([\d.]+)\.summary\.txt\.gz", listing)
    if not fn: raise LookupError(f"no summary file found under {sub}")
    ver = fn.group(1); files = {}
    for kind in ("summary", "changed_select_accessions", "protein_coding_genes_not_in_mane"):
        name = f"MANE.GRCh38.v{ver}.{kind}.txt.gz"; p = cache / name
        if not p.exists():
            try: urllib.request.urlretrieve(f"{FTP}/{sub}/{name}", p)
            except Exception: p = None
        files[kind] = p
    m = parse_summary(files["summary"], f"v{ver} (downloaded)"); _add_extras(m, files["changed_select_accessions"], files["protein_coding_genes_not_in_mane"]); return m
