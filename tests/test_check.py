import pathlib, pytest
from transcriptqc import load_mane, check_transcript, check_hgvs, check_file, mane_for_gene, summary, to_text, to_csv, to_json, to_html
from transcriptqc.mane import split_acc

M = load_mane()   # bundled v1.5


def test_snapshot_loaded():
    assert M.release.startswith("v1.5") and len(M.entries) > 19000 and sum(1 for e in M.entries if e.status == "MANE Plus Clinical") > 50
    assert split_acc("NM_007294.3") == ("NM_007294", 3) and split_acc("ENST00000357654") == ("ENST00000357654", None)


def test_gene_lookup_real_values():
    b = mane_for_gene("BRCA1", M); assert b["mane_select"] == "NM_007294.4" and b["mane_select_ensembl"].startswith("ENST00000357654")
    s = mane_for_gene("SCN2A", M); assert s["mane_select"] == "NM_001040142.2" and any(r.startswith("NM_001371246") for r, _ in s["plus_clinical"])
    assert not mane_for_gene("NOTAGENE", M)["found"]


def test_statuses():
    assert check_transcript("NM_007294.4", "BRCA1", M).status == "MANE_SELECT"
    r = check_transcript("NM_007294.3", "BRCA1", M); assert r.status == "OLD_VERSION" and r.mane_select == "NM_007294.4" and r.problem
    assert check_transcript("NM_007294.9", "BRCA1", M).status == "NEWER_VERSION"
    assert check_transcript("NM_007294", "BRCA1", M).status == "MANE_SELECT"
    r = check_transcript("NM_000059.4", "BRCA1", M); assert r.status == "GENE_MISMATCH" and r.gene_of_transcript == "BRCA2"
    r = check_transcript("NM_021007.2", "SCN2A", M); assert r.status == "NON_MANE" and "Plus Clinical" in r.note
    r = check_transcript("NM_001371246.1", "SCN2A", M); assert r.status == "MANE_PLUS_CLINICAL" and not r.problem
    assert check_transcript("NM_000000.1", "NOTAGENE", M).status == "UNKNOWN_GENE"
    assert check_transcript("NM_000000.1", None, M).status == "UNKNOWN_TRANSCRIPT"
    assert check_transcript("junk", "BRCA1", M).status == "UNPARSEABLE"
    assert check_transcript("ENST00000357654.9", "BRCA1", M).status == "MANE_SELECT" and check_transcript("ENST00000357654.5", None, M).status == "OLD_VERSION"


def test_changed_and_not_in_mane_from_real_release_files():
    r = check_transcript("NM_001302622.2", None, M); assert r.status == "MANE_SELECT_CHANGED" and r.gene_of_transcript == "ACHE" and "NM_000665.5" in r.note and r.mane_select == "NM_000665.5"
    assert "ACYP2" in {c["symbol"] for c in M.changed.values()} and M.changed[split_acc("NM_001320586.2")[0]]["affects_cds"] == "Yes"
    assert "AKR7L" in M.not_in_mane and check_transcript("NM_000000.1", "AKR7L", M).status == "GENE_NOT_IN_MANE"


def test_hgvs():
    r = check_hgvs("NM_007294.3:c.5266dupC", M); assert r.status == "OLD_VERSION" and r.variant == "c.5266dupC" and r.gene_of_transcript == "BRCA1"
    r = check_hgvs("NM_007294.3(BRCA1):c.5266dupC", M); assert r.gene == "BRCA1"
    assert check_hgvs("nonsense", M).status == "UNPARSEABLE"


def test_file_and_reports(tmp_path):
    res = check_file(pathlib.Path(__file__).parent.parent / "examples" / "variants.csv", M)
    st = {(r.gene, r.input): r.status for r in res}
    assert st[("BRCA1", "NM_007294.3")] == "OLD_VERSION" and st[("TP53", "NM_000546.6")] == "MANE_SELECT" and st[("SCN2A", "NM_021007.2")] == "NON_MANE" and st[("BRCA1", "NM_000059.4")] == "GENE_MISMATCH"
    s = summary(res); assert s["n"] == 14 and s["problems"] >= 5 and "SCN2A" in s["genes_with_plus_clinical"]
    assert "Potential reporting problems" in to_text(res, M.release) and to_csv(res).count("\n") == 15 and '"mane_release"' in to_json(res, M.release)
    h = to_html(res, M.release); assert "<table>" in h and "potential reporting problems" in h
    p = tmp_path / "h.tsv"; p.write_text("hgvs\nNM_000546.6:c.743G>A\n"); assert check_file(p, M)[0].status == "MANE_SELECT"
