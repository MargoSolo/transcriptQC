from __future__ import annotations
import csv, io, json, html, datetime as dt
from .check import Result, STATUS_HELP

ICON = {"MANE_SELECT": "✓", "MANE_PLUS_CLINICAL": "ℹ", "OLD_VERSION": "⚠", "OLD_VERSION_PLUS_CLINICAL": "⚠", "NEWER_VERSION": "⚠", "MANE_SELECT_CHANGED": "⚠", "NON_MANE": "⚠",
        "GENE_MISMATCH": "✗", "GENE_NOT_IN_MANE": "ℹ", "UNKNOWN_TRANSCRIPT": "?", "UNKNOWN_GENE": "?", "UNPARSEABLE": "✗"}
ORDER = ["GENE_MISMATCH", "UNPARSEABLE", "MANE_SELECT_CHANGED", "OLD_VERSION", "OLD_VERSION_PLUS_CLINICAL", "NON_MANE", "NEWER_VERSION", "UNKNOWN_TRANSCRIPT", "UNKNOWN_GENE", "GENE_NOT_IN_MANE", "MANE_PLUS_CLINICAL", "MANE_SELECT"]


def summary(results: list[Result]) -> dict:
    c = {s: 0 for s in ORDER}
    for r in results: c[r.status] = c.get(r.status, 0) + 1
    return {"n": len(results), "counts": {k: v for k, v in c.items() if v}, "problems": sum(1 for r in results if r.problem),
            "genes_with_plus_clinical": sorted({r.gene or r.gene_of_transcript for r in results if r.plus_clinical and (r.gene or r.gene_of_transcript)})}


def to_text(results: list[Result], release: str) -> str:
    L = []
    for r in results:
        head = f"{r.gene or r.gene_of_transcript or '?'}" + (f"  {r.variant}" if r.variant else "")
        L += [head, f"  Input transcript:      {r.input}", f"  Current MANE Select:   {r.mane_select or '—'}"]
        if r.plus_clinical: L.append(f"  MANE Plus Clinical:    {', '.join(r.plus_clinical)}")
        L += [f"  Status:                {r.status} {ICON.get(r.status, '')}", f"  Action:                {r.action}"]
        if r.note: L.append(f"  Note:                  {r.note}")
        L.append("")
    s = summary(results)
    L += ["TRANSCRIPT QC", "━" * 34, f"{s['n']} transcripts checked · MANE {release}", ""]
    for k, v in s["counts"].items(): L.append(f"{ICON[k]} {v:4d} {k.replace('_', ' ').lower()}")
    if s["genes_with_plus_clinical"]: L.append(f"ℹ {len(s['genes_with_plus_clinical']):4d} genes with a relevant MANE Plus Clinical transcript")
    L += ["", f"Potential reporting problems: {s['problems']}"]
    return "\n".join(L) + "\n"


def to_csv(results: list[Result]) -> str:
    b = io.StringIO(); w = csv.writer(b); w.writerow(["gene", "input", "variant", "mane_select", "mane_select_ensembl", "plus_clinical", "status", "action", "note"])
    for r in results: w.writerow([r.gene or r.gene_of_transcript or "", r.input, r.variant or "", r.mane_select or "", r.mane_select_ensembl or "", ";".join(r.plus_clinical), r.status, r.action, r.note])
    return b.getvalue()


def to_json(results: list[Result], release: str) -> str:
    return json.dumps({"mane_release": release, "summary": summary(results), "results": [r.to_dict() for r in results], "status_help": STATUS_HELP}, indent=2, ensure_ascii=False)


def to_html(results: list[Result], release: str, title="Transcript QC") -> str:
    s = summary(results); col = {"✓": "#2a9d8f", "⚠": "#e07a00", "✗": "#c0392b", "ℹ": "#1f5fa8", "?": "#7f8c8d"}
    css = "body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1100px;margin:30px auto;padding:0 20px;color:#222}h1{margin-bottom:2px}.sub{color:#666}.cards{display:flex;gap:12px;flex-wrap:wrap;margin:18px 0}.card{border:1px solid #e5e5e5;border-radius:10px;padding:12px 16px;min-width:150px}.card b{font-size:22px;display:block}table{border-collapse:collapse;font-size:13px;width:100%}th,td{border-bottom:1px solid #eee;padding:5px 8px;text-align:left;vertical-align:top}tr.p td{background:#fff8f0}code{font-family:Menlo,Consolas,monospace;font-size:12px}.foot{color:#888;font-size:12px;margin-top:30px}details{margin:14px 0}"
    H = [f"<title>{html.escape(title)}</title><style>{css}</style>", f"<h1>{html.escape(title)}</h1>", f"<p class=sub>{s['n']} transcripts checked against MANE {html.escape(release)} · {dt.date.today().isoformat()}</p>", "<div class=cards>"]
    for k, v in s["counts"].items(): H.append(f"<div class=card style='border-left:5px solid {col[ICON[k]]}'><b>{v}</b>{ICON[k]} {k.replace('_', ' ').lower()}</div>")
    H.append(f"<div class=card style='border-left:5px solid #c0392b'><b>{s['problems']}</b>potential reporting problems</div></div>")
    if s["genes_with_plus_clinical"]: H.append(f"<p>ℹ Genes with a MANE Plus Clinical transcript in this set: <b>{', '.join(html.escape(g) for g in s['genes_with_plus_clinical'])}</b> — check which isoform carries each variant.</p>")
    H += ["<table><tr><th></th><th>gene</th><th>input</th><th>variant</th><th>MANE Select</th><th>Plus Clinical</th><th>status</th><th>action</th><th>note</th></tr>"]
    for r in sorted(results, key=lambda r: ORDER.index(r.status)):
        H.append(f"<tr class='{'p' if r.problem else ''}'><td style='color:{col[ICON[r.status]]};font-weight:700'>{ICON[r.status]}</td><td>{html.escape(r.gene or r.gene_of_transcript or '')}</td><td><code>{html.escape(r.input)}</code></td><td><code>{html.escape(r.variant or '')}</code></td>"
                 f"<td><code>{html.escape(r.mane_select or '—')}</code></td><td><code>{html.escape(', '.join(r.plus_clinical))}</code></td><td>{r.status}</td><td>{html.escape(r.action)}</td><td>{html.escape(r.note)}</td></tr>")
    H.append("</table><details><summary>What the statuses mean</summary><ul>" + "".join(f"<li><b>{k}</b> — {html.escape(v)}</li>" for k, v in STATUS_HELP.items()) + "</ul></details>")
    H.append(f"<p class=foot>transcriptqc · MANE {html.escape(release)} from NCBI/EBI bulk files (offline, reproducible). A transcript check, not a variant classification.</p>")
    return "\n".join(H)
