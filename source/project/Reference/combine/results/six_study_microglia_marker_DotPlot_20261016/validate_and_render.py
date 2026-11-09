#!/usr/bin/python3
"""Server-side data checks and Poppler rendering; no source-object modifications."""
import csv, hashlib, json, math, pathlib, re, statistics, subprocess, sys
from PIL import Image, ImageOps, ImageDraw

root = pathlib.Path(sys.argv[1])
qa = root / "audit" / "pdf_render"
qa.mkdir(parents=True, exist_ok=True)
studies = ["Su", "Chia", "Fullard", "Jakel", "Grubman", "Yao"]
genes = "CSF1R P2RY12 AIF1 CX3CR1 TMEM119 OLFML3 GPR34 HEXB SALL1 GPNMB SPP1 APOE TREM2 CST3".split()
data_checks = []
render_checks = []
rendered = []
for study in studies:
    with (root / "tables" / f"{study}_microglia_markers_by_source_celltype.csv").open() as f:
        rows = list(csv.DictReader(f))
    labels = list(dict.fromkeys(row["celltype"] for row in rows))
    assert len(rows) == 14 * len(labels)
    assert {r["gene"] for r in rows} == set(genes)
    z_error = 0.0
    for gene in genes:
        group = [r for r in rows if r["gene"] == gene]
        assert len({r["celltype"] for r in group}) == len(labels)
        if group[0]["gene_present"] == "FALSE":
            assert all(r["percent_expressed"] == "NA" for r in group)
            continue
        values = [float(r["log1p_mean_CP10K"]) for r in group]
        mu = statistics.mean(values); sd = statistics.stdev(values)
        for row, value in zip(group, values):
            assert 0 <= float(row["percent_expressed"]) <= 100
            assert float(row["mean_normalized_CP10K"]) >= 0
            expected = min(2.5, max(-2.5, (value-mu)/sd))
            actual = float(row["scaled_average_expression"])
            assert math.isfinite(actual)
            z_error = max(z_error, abs(actual-expected))
    assert z_error < 1e-9
    data_checks.append({"study": study, "rows": len(rows), "celltypes": len(labels), "genes": 14, "max_independent_z_error": z_error})
    pdf = root / "pdf" / f"{study}_microglia_markers_source_celltype_DotPlot.pdf"
    info = subprocess.check_output(["/usr/bin/pdfinfo", str(pdf)], text=True)
    assert re.search(r"Pages:\s+1\b", info)
    # -raw preserves character order for 45-degree gene labels; -layout
    # interleaves characters from adjacent rotated words by their y position.
    text = subprocess.check_output(["/usr/bin/pdftotext", "-raw", str(pdf), "-"], text=True)
    flat = re.sub(r"\s+", "", text)
    assert all(g in flat for g in genes), (study, "missing gene text")
    assert all(re.sub(r"\s+", "", label) in flat for label in labels), (study, "missing celltype text")
    prefix = qa / study
    subprocess.run(["/usr/bin/pdftoppm", "-f", "1", "-singlefile", "-scale-to", "1800", "-png", str(pdf), str(prefix)], check=True)
    rendered.append(prefix.with_suffix(".png"))
    image = Image.open(prefix.with_suffix(".png"))
    assert image.width > 1000 and image.height > 400
    render_checks.append({"study":study,"pages":1,"all_14_gene_labels_found":True,"all_celltype_labels_found":True,"rendered_size":[image.width,image.height]})

combined = root / "Six_study_microglia_markers_source_celltype_DotPlots.pdf"
info = subprocess.check_output(["/usr/bin/pdfinfo", str(combined)], text=True)
assert re.search(r"Pages:\s+6\b", info)
subprocess.run(["/usr/bin/pdftoppm", "-scale-to", "1400", "-png", str(combined), str(qa / "Combined")], check=True)
all_pngs = sorted(qa.glob("Combined-*.png"))
assert len(all_pngs) == 6
thumbs = []
for p in rendered:
    im = Image.open(p).convert("RGB")
    im.thumbnail((950, 740))
    cell = Image.new("RGB", (970, 770), "#e8e8e8")
    cell.paste(im, ((970-im.width)//2, (770-im.height)//2))
    thumbs.append(cell)
contact = Image.new("RGB", (1940, 2310), "white")
for i, im in enumerate(thumbs): contact.paste(im, ((i%2)*970,(i//2)*770))
contact.save(qa / "All_six_single_PDF_contact_sheet.png")
thumbs = []
for p in all_pngs:
    im = Image.open(p).convert("RGB"); im.thumbnail((900,680))
    cell = Image.new("RGB",(920,700),"#e8e8e8"); cell.paste(im,((920-im.width)//2,(700-im.height)//2)); thumbs.append(cell)
contact = Image.new("RGB",(1840,2100),"white")
for i, im in enumerate(thumbs): contact.paste(im,((i%2)*920,(i//2)*700))
contact.save(qa / "Combined_PDF_contact_sheet.png")
assert len(list((root/"pdf").glob("*.pdf")))==6
assert len(list((root/"png").glob("*.png")))==6
checks = {"data_checks":data_checks,"individual_pdf_checks":render_checks,"combined_pdf_pages":6,"independent_zscore_validation_passed":True,"visual_review_status":"pending manual visual inspection of rendered PDFs"}
(root / "audit" / "validation_report.json").write_text(json.dumps(checks,indent=2))
hashes = {}
for p in sorted(list((root/"pdf").glob("*.pdf"))+list((root/"png").glob("*.png"))+[combined]):
    hashes[str(p.relative_to(root))] = hashlib.sha256(p.read_bytes()).hexdigest()
(root/"audit"/"final_figure_sha256.json").write_text(json.dumps(hashes,indent=2))
print(json.dumps(checks,indent=2))
