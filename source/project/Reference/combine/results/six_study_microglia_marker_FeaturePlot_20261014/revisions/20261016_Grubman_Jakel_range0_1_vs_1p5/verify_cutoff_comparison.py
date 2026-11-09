#!/usr/bin/python3
"""Validate only the Grubman/Jakel cutoff comparison and protected Fullard files.

Run on server 24:
  /usr/bin/python3 verify_cutoff_comparison.py [output_dir] [--baseline v1_dir]
No expression analysis is performed. All generated files remain in output_dir.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone

from PIL import Image


STUDIES = ("Grubman", "Jakel")
MARKERS = ("P2RY12", "CX3CR1", "CSF1R")
VERSIONS = {"range0_1p5": 1.5, "range0_1": 1.0}
BASELINE_NAME = "20261016_Grubman_Jakel_Fullard_range0_2"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(row, key):
    return float(row[key])


def integer(row, key):
    value = number(row, key)
    if not math.isfinite(value) or not value.is_integer():
        raise ValueError("Invalid integer {}={!r}".format(key, row[key]))
    return int(value)


def close(left, right):
    return math.isclose(float(left), float(right), rel_tol=1e-9, abs_tol=1e-8)


def truth(value):
    return str(value).strip().lower() in ("true", "t", "1")


def run_command(arguments):
    result = subprocess.run(arguments, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=False)
    if result.returncode:
        raise RuntimeError("{} failed ({}): {}".format(
            arguments[0], result.returncode, result.stderr.strip()))
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", nargs="?", default=str(Path(__file__).resolve().parent))
    parser.add_argument("--baseline")
    args = parser.parse_args()
    root = Path(args.output_dir).expanduser().resolve()
    baseline = (Path(args.baseline).expanduser().resolve() if args.baseline
                else root.parent / BASELINE_NAME)
    audit = root / "audit"
    audit.mkdir(parents=True, exist_ok=True)
    report = {
        "status": "running", "revision": str(root), "baseline": str(baseline),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "validation_source": "Current files, baseline audits, protected hash snapshot and fresh Poppler renders",
        "errors": [],
    }

    def check(condition, message):
        if not condition:
            report["errors"].append(message)
        return bool(condition)

    def section(name, operation):
        try:
            operation()
        except Exception as exc:
            report["errors"].append("{}: {}: {}".format(name, type(exc).__name__, exc))

    expected = {}
    for version in VERSIONS:
        for study in STUDIES:
            for marker in MARKERS:
                stem = "{}_{}_FeaturePlot".format(study, marker)
                expected[version + "/pdf_single_marker/" + stem + ".pdf"] = {
                    "labels": [study, marker], "inches": (4.8, 4.6),
                    "png": version + "/png_single_marker/" + stem + ".png",
                    "png_pixels": (1440, 1380),
                }
            stem = study + "_three_microglia_markers_FeaturePlots"
            expected[version + "/pdf_by_study/" + stem + ".pdf"] = {
                "labels": [study] + list(MARKERS), "inches": (14.4, 5.0),
                "png": version + "/png_by_study/" + stem + ".png",
                "png_pixels": (3456, 1200),
            }
    for study in STUDIES:
        stem = study + "_range0_1p5_vs_1_FeaturePlots"
        expected["comparison/" + stem + ".pdf"] = {
            "labels": [study] + list(MARKERS), "inches": (14.4, 9.6),
            "png": "comparison/" + stem + ".png", "png_pixels": (3456, 2304),
        }

    def verify_inventory():
        pdfs, pngs = set(), set()
        for version in VERSIONS:
            for folder in ("pdf_single_marker", "pdf_by_study"):
                pdfs.update(str(p.relative_to(root)) for p in (root / version / folder).glob("*.pdf"))
            for folder in ("png_single_marker", "png_by_study"):
                pngs.update(str(p.relative_to(root)) for p in (root / version / folder).glob("*.png"))
        pdfs.update(str(p.relative_to(root)) for p in (root / "comparison").glob("*.pdf"))
        pngs.update(str(p.relative_to(root)) for p in (root / "comparison").glob("*.png"))
        wanted_pngs = {item["png"] for item in expected.values()}
        check(pdfs == set(expected), "PDF inventory mismatch; missing={}, extra={}".format(
            sorted(set(expected) - pdfs), sorted(pdfs - set(expected))))
        check(pngs == wanted_pngs, "PNG inventory mismatch; missing={}, extra={}".format(
            sorted(wanted_pngs - pngs), sorted(pngs - wanted_pngs)))
        report["inventory"] = {"n_pdfs": len(pdfs), "n_pngs": len(pngs),
                               "expected_each": 18, "png_checks": []}
        for item in expected.values():
            path = root / item["png"]
            if not path.is_file():
                continue
            with Image.open(path) as img:
                size = img.size
                img.verify()
            check(size == item["png_pixels"], "PNG dimensions changed: {} has {} expected {}".format(
                item["png"], size, item["png_pixels"]))
            report["inventory"]["png_checks"].append({
                "path": item["png"], "size_pixels": size, "readable": True,
                "sha256": digest(path),
            })

    def verify_scales():
        rows = read_tsv(audit / "plot_scale_parameters.tsv")
        old_rows = read_tsv(baseline / "audit/plot_scale_parameters.tsv")
        old = {(row["study"], row["marker"]): row for row in old_rows}
        wanted = {(study, marker, version) for study in STUDIES for marker in MARKERS for version in VERSIONS}
        keyed = {}
        for row in rows:
            key = (row["study"], row["marker"], row["version"])
            check(key not in keyed, "Duplicate scale record: {}".format(key))
            keyed[key] = row
        check(len(rows) == 12 and set(keyed) == wanted, "Expected exactly 12 target scale records")
        checked = []
        for key in sorted(wanted & set(keyed)):
            study, marker, version = key
            row = keyed[key]
            previous = old[(study, marker)]
            upper = VERSIONS[version]
            prefix = "{}/{}/{}".format(study, marker, version)
            check(close(number(row, "lower_limit"), 0), prefix + ": lower limit is not 0")
            check(close(number(row, "upper_limit"), upper), prefix + ": wrong upper limit")
            ticks = [float(value) for value in row["legend_breaks"].split(";")]
            expected_ticks = [0, 0.5, 1, 1.5] if upper == 1.5 else [0, 0.5, 1]
            check(len(ticks) == len(expected_ticks) and all(close(a, b) for a, b in zip(ticks, expected_ticks)),
                  prefix + ": wrong legend ticks")
            check(row["colour_low"].upper() == "#D3D3D3" and row["colour_high"].upper() == "#E31A1C",
                  prefix + ": plotting colours changed")
            for name, value in (("pt_size", 1.3), ("raster_resolution_x", 1200),
                                ("raster_resolution_y", 1200), ("coordinate_ratio", 1)):
                check(close(number(row, name), value), prefix + ": " + name + " changed")
            for name in ("n_cells", "n_positive"):
                check(integer(row, name) == integer(previous, name), prefix + ": " + name + " changed")
            n_positive, n_saturated = integer(row, "n_positive"), integer(row, "n_saturated")
            check(0 <= n_saturated <= n_positive <= integer(row, "n_cells"), prefix + ": invalid counts")
            pct = 100 * n_saturated / max(1, n_positive)
            check(close(number(row, "percent_saturated_of_positive"), pct), prefix + ": saturation percentage mismatch")
            checked.append(dict(row))
        monotonic = []
        for study in STUDIES:
            for marker in MARKERS:
                pair = (study, marker)
                if all(pair + (version,) in keyed for version in VERSIONS):
                    one = integer(keyed[pair + ("range0_1",)], "n_saturated")
                    one_half = integer(keyed[pair + ("range0_1p5",)], "n_saturated")
                    two = integer(old[pair], "n_saturated")
                    check(close(number(old[pair], "upper_limit"), 2), str(pair) + ": baseline upper limit is not 2")
                    check(one >= one_half >= two, str(pair) + ": saturation counts not monotonic: 1 >= 1.5 >= 2")
                    monotonic.append({"study": study, "marker": marker, "at_1": one,
                                      "at_1p5": one_half, "at_2": two})
        report["scales"] = {"n_records": len(rows), "records": checked,
                            "saturation_monotonicity": monotonic}

    def verify_coordinates():
        rows = read_tsv(audit / "plot_coordinate_validation.tsv")
        check(len(rows) == 2 and {row["study"] for row in rows} == set(STUDIES),
              "Coordinate audit must have exactly Grubman and Jakel")
        old = {(row["study"], row["marker"]): row
               for row in read_tsv(baseline / "audit/plot_scale_parameters.tsv")}
        for row in rows:
            study = row["study"]
            check(integer(row, "n_cells") == integer(old[(study, MARKERS[0])], "n_cells"),
                  study + ": coordinate audit cell count changed")
            for name in ("coordinate_values_identical", "barcode_order_identical",
                         "expression_values_identical_to_plot_bundle"):
                check(truth(row[name]), study + ": " + name + " is not TRUE")
            check(close(number(row, "coordinate_ratio"), 1), study + ": coordinate ratio changed")
            check(str(row["umap_recomputed"]).strip().lower() in ("false", "f", "0"),
                  study + ": UMAP recomputation is not FALSE")
        report["coordinates"] = rows

    def verify_input_hashes():
        checks = []
        for study in STUDIES:
            rel = "plot_data/" + study + "_stored_umap_three_RNA_markers.rds"
            current_hash, old_hash = digest(root / rel), digest(baseline / rel)
            same = current_hash == old_hash
            check(same, "Plot bundle changed: " + rel)
            checks.append({"path": rel, "sha256": current_hash,
                           "baseline_sha256": old_hash, "identical": same})
        report["plot_data_hashes"] = checks

    def verify_fullard():
        snapshot_path = audit / "protected_Fullard_before.json"
        payload = json.loads(snapshot_path.read_text(encoding="utf-8"))
        check(Path(payload["approved_revision"]).resolve() == baseline.resolve(),
              "Protected Fullard snapshot refers to a different baseline")
        files = payload["fullard_files"]
        check(isinstance(files, dict) and len(files) == 12,
              "Protected Fullard snapshot must contain all 12 approved PNG/PDF files")
        checks = []
        for relative, before_hash in sorted(files.items()):
            rel = Path(relative)
            check(not rel.is_absolute() and ".." not in rel.parts,
                  "Protected file must be relative to baseline: " + relative)
            check("Fullard" in rel.name and rel.suffix.lower() in (".png", ".pdf"),
                  "Unexpected protected Fullard file: " + relative)
            path = baseline / rel
            after_hash = digest(path)
            same = after_hash == before_hash
            check(same, "Approved Fullard file changed: " + relative)
            checks.append({"path": relative, "before_sha256": before_hash,
                           "after_sha256": after_hash, "identical": same})
        scale_rows = read_tsv(baseline / "audit/plot_scale_parameters.tsv")
        fullard = [row for row in scale_rows if row["study"] == "Fullard"]
        check(len(fullard) == 3 and {row["marker"] for row in fullard} == set(MARKERS),
              "Approved Fullard scale audit must contain the three markers")
        for row in fullard:
            check(close(number(row, "lower_limit"), 0) and close(number(row, "upper_limit"), 2),
                  "Approved Fullard scale changed: " + row["marker"])
        report["protected_Fullard"] = {"n_files": len(files), "files": checks,
                                      "approved_scale_records": fullard}

    def verify_pdfs():
        tools = {name: shutil.which(name) for name in ("pdfinfo", "pdftotext", "pdftoppm")}
        report["tools"] = tools
        if not all(check(value is not None, "Required Poppler tool missing: " + name)
                   for name, value in tools.items()):
            return
        checks = []
        for relative, item in sorted(expected.items()):
            path = root / relative
            if not check(path.is_file(), "Expected PDF missing: " + relative):
                continue
            info = run_command([tools["pdfinfo"], str(path)])
            match = re.search(r"^Pages:\s+(\d+)\s*$", info, flags=re.MULTILINE)
            pages = int(match.group(1)) if match else None
            check(pages == 1, "PDF must contain exactly one page: " + relative)
            text = run_command([tools["pdftotext"], "-layout", str(path), "-"])
            for label in item["labels"]:
                check(label in text, "PDF title label missing: {} in {}".format(label, relative))
            page_size = re.search(r"^Page size:\s+([\d.]+)\s+x\s+([\d.]+)\s+pts", info, flags=re.MULTILINE)
            if check(page_size is not None, "Cannot inspect PDF page size: " + relative):
                pts = tuple(float(value) for value in page_size.groups())
                # cairo_pdf quantizes page dimensions to whole PostScript points.
                check(all(abs(actual - inch * 72) < 1 for actual, inch in zip(pts, item["inches"])),
                      "PDF page size changed: " + relative)
                if not relative.startswith("comparison/"):
                    baseline_pdf = baseline / Path(relative).relative_to(Path(relative).parts[0])
                    old_info = run_command([tools["pdfinfo"], str(baseline_pdf)])
                    old_size = re.search(r"^Page size:\s+([\d.]+)\s+x\s+([\d.]+)\s+pts", old_info, flags=re.MULTILINE)
                    check(old_size is not None and pts == tuple(float(v) for v in old_size.groups()),
                          "PDF dimensions differ from approved baseline: " + relative)
            render_prefix = root / "rendered_pdf_check" / Path(relative).with_suffix("")
            render_prefix.parent.mkdir(parents=True, exist_ok=True)
            run_command([tools["pdftoppm"], "-r", "100", "-png", "-singlefile", str(path), str(render_prefix)])
            render_path = render_prefix.with_suffix(".png")
            with Image.open(render_path) as img:
                rendered_size = img.size
                img.verify()
            check(all(abs(actual - inch * 100) <= 1 for actual, inch in zip(rendered_size, item["inches"])),
                  "Unexpected rendered PDF size: " + relative)
            if relative.startswith("comparison/"):
                check("Top row: 0-1.5; bottom row: 0-1." in text,
                      "Comparison row/cutoff explanation missing: " + relative)
            checks.append({"pdf": relative, "pages": pages, "labels_verified": item["labels"],
                           "render_dpi": 100, "rendered_path": str(render_path.relative_to(root)),
                           "size_pixels": rendered_size, "sha256": digest(path)})
        report["pdf_validation"] = {"n_pdfs": len(checks), "n_total_pages": sum(
            item["pages"] or 0 for item in checks), "render_dpi": 100, "pdfs": checks}
        check(len(checks) == 18, "Fresh PDF validation/rendering did not complete for all 18 PDFs")

    for name, operation in (("inventory", verify_inventory), ("scales", verify_scales),
                            ("coordinates", verify_coordinates), ("plot_data_hashes", verify_input_hashes),
                            ("protected_Fullard", verify_fullard), ("pdf_validation", verify_pdfs)):
        section(name, operation)
    report["status"] = "passed" if not report["errors"] else "failed"
    output = audit / "comparison_validation.json"
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("COMPARISON_VALIDATION_{}: {}".format(report["status"].upper(), output))
    for error in report["errors"]:
        print("ERROR: " + error, file=sys.stderr)
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
