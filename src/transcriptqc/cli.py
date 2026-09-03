from __future__ import annotations
import sys
from pathlib import Path
import typer
from .mane import load_mane
from .check import check_transcript, check_hgvs, check_file, mane_for_gene
from .report import to_text, to_csv, to_json, to_html, summary

app = typer.Typer(add_completion=False)


@app.command(help="Is this transcript right for clinical reporting? TARGET: a CSV/TSV file, a transcript (NM_…/ENST…), an HGVS string (NM_…:c.…) or a gene symbol.")
def main(target: str = typer.Argument(..., help="file · transcript · HGVS · gene"),
         mane: str = typer.Option("bundled", help="'bundled' (v1.5 snapshot), a release like '1.5', or 'current' (downloads once)"),
         csv_out: Path | None = typer.Option(None, "--csv"), json_out: Path | None = typer.Option(None, "--json"), html_out: Path | None = typer.Option(None, "--html"),
         gene: str | None = typer.Option(None, help="gene symbol when checking a single transcript"), fail: bool = typer.Option(False, help="exit 1 if any potential problem")):
    m = load_mane(mane); p = Path(target)
    if p.exists(): res = check_file(p, m)
    elif ":" in target: res = [check_hgvs(target, m)]
    elif target.upper().startswith(("NM_", "NR_", "XM_", "XR_", "ENST")): res = [check_transcript(target, gene, m)]
    else:
        info = mane_for_gene(target, m)
        if not info["found"] and not info["not_in_mane"]: typer.echo(f"{target}: not found in MANE {m.release}"); raise typer.Exit(2)
        typer.echo(f"{info['gene']} · MANE {info['release']}\n  MANE Select:        {info['mane_select']}  ({info['mane_select_ensembl']}; {info['mane_select_protein']})")
        for r, e in info["plus_clinical"]: typer.echo(f"  MANE Plus Clinical: {r}  ({e})")
        if info["not_in_mane"]: typer.echo("  protein-coding gene not yet in MANE")
        raise typer.Exit()
    typer.echo(to_text(res, m.release))
    if csv_out: csv_out.write_text(to_csv(res), encoding="utf-8")
    if json_out: json_out.write_text(to_json(res, m.release), encoding="utf-8")
    if html_out: html_out.write_text(to_html(res, m.release), encoding="utf-8")
    if fail and summary(res)["problems"]: sys.exit(1)


if __name__ == "__main__":
    app()
