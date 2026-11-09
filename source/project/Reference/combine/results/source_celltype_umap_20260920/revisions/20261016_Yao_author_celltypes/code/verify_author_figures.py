#!/usr/bin/env python3
"""Verify annotation audits, rendered figures and unchanged non-Yao pages."""
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from PIL import Image

task = Path(sys.argv[1]).resolve()
baseline = Path(__import__("os").path.join(__import__("os").environ["MG_PROJECT_ROOT"], "Reference/combine/results/source_celltype_umap_20260920/revisions/20261016_celltype_labels_on_umap_no_boxes"))
checks = []
def check(name, condition, details=None):
    checks.append({"check": name, "passed": bool(condition), "details": details})
    if not condition:
        raise AssertionError(name)
def command(*args):
    return subprocess.check_output([str(v) for v in args], text=True)
def raster(pdf, folder):
    folder.mkdir(parents=True, exist_ok=True)
    subprocess.run(["pdftoppm","-r","100","-png",str(pdf),str(folder/"page")],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    return sorted(folder.glob("page-*.png"))
def pixel_hash(path):
    with Image.open(path) as im:
        return (im.size, hashlib.sha256(im.convert("RGB").tobytes()).hexdigest())

mapping = json.loads((task/"audit/annotation_mapping_validation.json").read_text())
check("73721 retained, 73547 author labels, 174 explicitly unmatched",
      (mapping["target_n"],mapping["matched_n"],mapping["unmatched_n"]) == (73721,73547,174))
check("unique sample+barcode mapping and independent author class counts",
      mapping["sample_and_barcode_keys_unique"] and mapping["original_16_class_counts_match_paper"] and mapping["current_16_class_counts_independently_verified"])
check("saved coordinates preserved without fitting", mapping["current_coordinates_identical"] and mapping["original_coordinates_identical"] and not mapping["umap_recomputed"])
check("unmatched labels were not inferred", not mapping["unmatched_labels_inferred"])
check("cell identity corroborated by original count metrics", mapping["matched_QC_both_exact"] == 60044)
check("author source did not change during extraction",
      json.loads((task/"audit/author_object_inventory.json").read_text())["source_unchanged_size_mtime"])
counts = list(csv.DictReader((task/"audit/plot_counts_and_colours.csv").open()))
specs = [("Yao_current_UMAP_author_labels",1),("Yao_author_atlas_original_UMAP_labels",1),("Six_study_celltypes_with_Yao_author_labels",6)]
rendered = {}
for plot_id, expected_pages in specs:
    pdf = task/"output/pdf"/(plot_id+".pdf")
    info = command("pdfinfo",pdf)
    pages = int(re.search(r"^Pages:\s+(\d+)",info,re.M).group(1))
    check(plot_id+": PDF page count", pages==expected_pages)
    check(plot_id+": PDF dimensions",
          bool(re.search(r"Page size:\s+900(?:\.\d+)? x 720(?:\.\d+)? pts",info)))
    rendered[plot_id] = raster(pdf,task/"audit/pdf_render"/plot_id)
    check(plot_id+": every PDF page rendered",len(rendered[plot_id])==expected_pages)
    if expected_pages==1:
        text = command("pdftotext","-layout",pdf,"-")+command("pdftotext","-raw",pdf,"-")
        compact = re.sub(r"\s+","",text)
        for row in (x for x in counts if x["plot_id"]==plot_id):
            check(plot_id+": label "+row["display_label"], re.sub(r"\s+","",row["display_label"]) in compact)
            check(plot_id+": count "+row["n_cells"], "("+format(int(row["n_cells"]),",")+")" in compact)
        png = task/"output/png"/(plot_id+".png")
        with Image.open(png) as im:
            check(plot_id+": PNG resolution",im.size==(3000,2400) and all(abs(v-240)<0.1 for v in im.info["dpi"]))
old = raster(baseline/"Six_study_celltypes_labeled_on_UMAPs.pdf",task/"audit/pdf_render/baseline_main")
new = rendered["Six_study_celltypes_with_Yao_author_labels"]
for i in range(5):
    check("non-Yao page "+str(i+1)+": identical rendered pixels",pixel_hash(old[i])==pixel_hash(new[i]))
check("six-study book Yao page equals updated independent figure",
      pixel_hash(new[5])==pixel_hash(rendered["Yao_current_UMAP_author_labels"][0]))
result = {"passed":all(x["passed"] for x in checks),"n_checks":len(checks),"n_pdfs":3,"n_pdf_pages":8,"checks":checks}
(task/"audit/figure_validation.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
print(json.dumps({k:v for k,v in result.items() if k!="checks"}))
