#!/usr/bin/python3
"""Verify this display-only revision against the unmodified plotting bundle.

Usage: /usr/bin/python3 verify_revision.py [revision_dir] [--baseline baseline_dir]
Core checks use only Python's standard library and Poppler command-line tools.
Pillow is optional and is used only to assemble the nine-panel contact sheet.
No packages are installed and no source data or figures are modified.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
from datetime import datetime, timezone


STUDIES = ("Su", "Chia", "Fullard", "Jakel", "Grubman", "Yao")
MARKERS = ("P2RY12", "CX3CR1", "CSF1R")
TARGET_STUDIES = ("Grubman", "Jakel", "Fullard")
AUDIT_BOOK = "Six_study_existing_MG_annotation_and_marker_audit.pdf"
MAIN_FIGURE = "Six_study_P2RY12_CX3CR1_CSF1R_FeaturePlots.pdf"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha256(path):
    require(path.is_file(), f"Required file is missing: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_scale_rows(root):
    path = root / "audit" / "plot_scale_parameters.tsv"
    require(path.is_file(), f"Scale audit is missing: {path}")
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    require(len(rows) == 18, f"Expected 18 scale records in {path}, found {len(rows)}")
    expected = {(study, marker) for study in STUDIES for marker in MARKERS}
    keyed = {}
    for row in rows:
        key = (row["study"], row["marker"])
        require(key not in keyed, f"Duplicate scale record in {path}: {key}")
        keyed[key] = row
    require(set(keyed) == expected, f"Unexpected study/marker records in {path}")
    return keyed


def number_equal(left, right):
    left, right = float(left), float(right)
    return math.isfinite(left) and math.isfinite(right) and math.isclose(
        left, right, rel_tol=1e-12, abs_tol=1e-12
    )


def verify_scales(revision, baseline, report):
    revised = read_scale_rows(revision)
    original = read_scale_rows(baseline)
    results = []
    changed_keys = set()
    for study in STUDIES:
        for marker in MARKERS:
            key = (study, marker)
            row, old = revised[key], original[key]
            target = study in TARGET_STUDIES
            low, high = float(row["lower_limit"]), float(row["upper_limit"])
            if target:
                require(number_equal(low, 0) and number_equal(high, 2),
                        f"Target scale must be 0-2: {study}/{marker}, got {low}-{high}")
            else:
                for field in ("lower_limit", "upper_limit"):
                    require(number_equal(row[field], old[field]),
                            f"Non-target {study}/{marker} changed {field}: {old[field]} -> {row[field]}")
            if any(not number_equal(row[field], old[field])
                   for field in ("lower_limit", "upper_limit")):
                changed_keys.add(key)
            # These rendering and biological counts must be unchanged for every panel.
            for field in ("colour_low", "colour_high", "n_cells", "n_positive",
                          "pt_size", "raster_resolution_x", "raster_resolution_y",
                          "coordinate_ratio"):
                require(row[field] == old[field],
                        f"Unexpected change in {study}/{marker} field {field}: {old[field]} -> {row[field]}")
            if not target:
                for field in ("n_saturated", "percent_saturated_of_positive"):
                    require(number_equal(row[field], old[field]),
                            f"Non-target saturation changed: {study}/{marker}/{field}")
            results.append({"study": study, "marker": marker, "target": target,
                            "lower_limit": low, "upper_limit": high,
                            "baseline_lower_limit": float(old["lower_limit"]),
                            "baseline_upper_limit": float(old["upper_limit"])})
    targets = {(study, marker) for study in TARGET_STUDIES for marker in MARKERS}
    require(changed_keys == targets,
            f"Exactly the nine target scales must change; changed keys: {sorted(changed_keys)}")
    report["scales"] = {"n_records": 18, "n_target_0_2": 9,
                        "n_non_target_unchanged": 9, "records": results}


def verify_hashes(revision, baseline, report):
    data_results = []
    expected_data = {f"{study}_stored_umap_three_RNA_markers.rds" for study in STUDIES}
    for root in (revision, baseline):
        found = {path.name for path in (root / "plot_data").glob("*.rds")}
        require(found == expected_data, f"Expected exactly six plot_data RDS files in {root}; found {sorted(found)}")
    for name in sorted(expected_data):
        relative = Path("plot_data") / name
        new_hash, old_hash = sha256(revision / relative), sha256(baseline / relative)
        require(new_hash == old_hash, f"Plot data changed: {relative}")
        data_results.append({"path": str(relative), "sha256": new_hash,
                             "baseline_sha256": old_hash, "identical": True})
    report["plot_data_hashes"] = data_results
    png_results = []
    for study in STUDIES:
        if study in TARGET_STUDIES:
            continue
        for marker in MARKERS:
            relative = Path("png_single_marker") / f"{study}_{marker}_FeaturePlot.png"
            new_hash, old_hash = sha256(revision / relative), sha256(baseline / relative)
            bytes_identical = new_hash == old_hash
            pixels_identical = None
            pixel_check = "not_required_bytes_identical"
            if not bytes_identical:
                try:
                    from PIL import Image
                except ImportError as error:
                    raise AssertionError(
                        f"Non-target PNG byte hashes differ and Pillow is unavailable for pixel comparison: {relative}"
                    ) from error
                with Image.open(revision / relative) as new_image, Image.open(baseline / relative) as old_image:
                    new_rgb, old_rgb = new_image.convert("RGBA"), old_image.convert("RGBA")
                    pixels_identical = (new_rgb.size == old_rgb.size and new_rgb.tobytes() == old_rgb.tobytes())
                pixel_check = "performed_after_byte_hash_mismatch"
            png_results.append({"path": str(relative), "sha256": new_hash,
                                "baseline_sha256": old_hash, "bytes_identical": bytes_identical,
                                "pixels_identical": pixels_identical, "pixel_check": pixel_check,
                                "visually_identical": bytes_identical or pixels_identical is True})
            report["non_target_png_hashes"] = png_results
            require(bytes_identical or pixels_identical is True,
                    f"Non-target PNG pixels changed: {relative}")
    require(len(png_results) == 9, "Expected nine unchanged non-target single-marker PNGs")
    report["non_target_png_hashes"] = png_results


def run_tool(arguments):
    result = subprocess.run(arguments, check=False, capture_output=True, text=True,
                            encoding="utf-8", errors="replace", timeout=180)
    require(result.returncode == 0,
            f"Command failed ({result.returncode}): {arguments!r}\n{result.stderr.strip()}")
    return result.stdout


def png_size(path):
    require(path.is_file(), f"Rendered PNG is missing: {path}")
    with path.open("rb") as stream:
        header = stream.read(24)
    require(len(header) == 24 and header[:8] == b"\x89PNG\r\n\x1a\n" and header[12:16] == b"IHDR",
            f"Invalid PNG header: {path}")
    width, height = struct.unpack(">II", header[16:24])
    require(width > 350 and height > 350, f"Unexpectedly small rendered PDF page: {path}, {width}x{height}")
    return [width, height]


def verify_and_render_pdfs(revision, report):
    tools = {name: shutil.which(name) for name in ("pdfinfo", "pdftotext", "pdftoppm")}
    require(all(tools.values()), f"Required Poppler tools are unavailable: {tools}")
    report["tools"] = tools
    expected = {Path(MAIN_FIGURE), Path(AUDIT_BOOK)}
    expected.update(Path("pdf_by_study") / f"{study}_three_microglia_markers_FeaturePlots.pdf"
                    for study in STUDIES)
    expected.update(Path("pdf_single_marker") / f"{study}_{marker}_FeaturePlot.pdf"
                    for study in STUDIES for marker in MARKERS)
    found = {path.relative_to(revision) for path in revision.glob("*.pdf")}
    found.update(path.relative_to(revision) for path in (revision / "pdf_by_study").glob("*.pdf"))
    found.update(path.relative_to(revision) for path in (revision / "pdf_single_marker").glob("*.pdf"))
    require(found == expected,
            f"PDF inventory mismatch; expected 26. Missing: {sorted(map(str, expected - found))}; "
            f"unexpected: {sorted(map(str, found - expected))}")
    render_dir = revision / "rendered_pdf_check"
    render_dir.mkdir(exist_ok=True)
    results = []
    for relative in sorted(expected, key=str):
        pdf = revision / relative
        info = run_tool([tools["pdfinfo"], str(pdf)])
        match = re.search(r"^Pages:\s*(\d+)\s*$", info, re.MULTILINE)
        require(match is not None, f"Could not determine page count: {relative}")
        pages = int(match.group(1))
        expected_pages = 6 if pdf.name == AUDIT_BOOK else 1
        require(pages == expected_pages,
                f"Wrong page count: {relative}, expected {expected_pages}, found {pages}")
        text = run_tool([tools["pdftotext"], "-layout", str(pdf), "-"])
        if relative.parent.name == "pdf_single_marker":
            study, marker, _ = pdf.stem.split("_", 2)
            required = [study, marker]
        elif relative.parent.name == "pdf_by_study":
            required = [pdf.name.split("_", 1)[0], *MARKERS]
        else:
            required = [*STUDIES, *MARKERS]
        require(all(label in text for label in required),
                f"Missing study/marker text in {relative}: {[label for label in required if label not in text]}")
        if pdf.name == AUDIT_BOOK:
            text_pages = [page for page in text.split("\f") if page.strip()]
            require(len(text_pages) == 6, f"Audit book text did not contain six pages: {relative}")
            for study, page_text in zip(STUDIES, text_pages):
                require(all(label in page_text for label in (study, *MARKERS)),
                        f"Missing study/marker label on {study} audit page")
        prefix = render_dir / pdf.stem
        # Remove only this PDF's prior render pages, so stale pages cannot pass checks.
        for old_page in render_dir.glob(pdf.stem + "-*.png"):
            if re.fullmatch(re.escape(pdf.stem) + r"-\d+\.png", old_page.name):
                old_page.unlink()
        run_tool([tools["pdftoppm"], "-r", "100", "-png", str(pdf), str(prefix)])
        rendered = sorted(render_dir.glob(pdf.stem + "-*.png"))
        require(len(rendered) == expected_pages,
                f"Rendered page count mismatch: {relative}, found {len(rendered)}")
        results.append({"pdf": str(relative), "pages": pages,
                        "labels_verified": required, "render_dpi": 100,
                        "rendered_pages": [{"path": str(path.relative_to(revision)),
                                            "size_pixels": png_size(path)} for path in rendered]})
    total_pages = sum(item["pages"] for item in results)
    require(len(results) == 26 and total_pages == 31,
            f"Expected 26 PDFs / 31 pages, found {len(results)} PDFs / {total_pages} pages")
    report["pdf_validation"] = {"n_pdfs": len(results), "n_total_pages": total_pages,
                                "render_dpi": 100, "pdfs": results}


def make_contact_sheet(revision, report):
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        report["contact_sheet"] = {"status": "skipped", "reason": "Pillow is unavailable; no packages installed"}
        return
    render_dir = revision / "rendered_pdf_check"
    width, height, title_height = 480, 460, 26
    canvas = Image.new("RGB", (3 * width, 3 * (height + title_height)), "white")
    draw = ImageDraw.Draw(canvas)
    sources = []
    for row, study in enumerate(TARGET_STUDIES):
        for col, marker in enumerate(MARKERS):
            page = render_dir / f"{study}_{marker}_FeaturePlot-1.png"
            with Image.open(page) as source:
                image = source.convert("RGB")
                image.thumbnail((width, height))
                x = col * width + (width - image.width) // 2
                y = row * (height + title_height) + title_height + (height - image.height) // 2
                canvas.paste(image, (x, y))
            draw.text((col * width + 10, row * (height + title_height) + 7),
                      f"{study} | {marker} | colour range 0-2", fill="black")
            sources.append(str(page.relative_to(revision)))
    path = render_dir / "Grubman_Jakel_Fullard_9_target_panels_contact_sheet.png"
    canvas.save(path)
    report["contact_sheet"] = {"status": "created", "path": str(path.relative_to(revision)),
                               "size_pixels": list(canvas.size), "sources": sources}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("revision", nargs="?", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--baseline", type=Path)
    arguments = parser.parse_args()
    revision = arguments.revision.resolve()
    baseline = (arguments.baseline or revision.parent.parent).resolve()
    report_path = revision / "audit" / "revision_validation.json"
    report = {"status": "running", "revision": str(revision), "baseline": str(baseline),
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "validation_source": "Current revision files and fresh Poppler renders; prior validation JSON is not reused",
              "prior_pdf_validation_json_present": (revision / "audit" / "local_pdf_validation.json").is_file()}
    try:
        require(revision != baseline, "Revision and baseline must be different directories")
        require(revision.is_dir(), f"Revision directory is missing: {revision}")
        require(baseline.is_dir(), f"Baseline directory is missing: {baseline}")
        verify_scales(revision, baseline, report)
        verify_hashes(revision, baseline, report)
        verify_and_render_pdfs(revision, report)
        make_contact_sheet(revision, report)
        report["status"] = "passed"
    except Exception as error:
        report["status"] = "failed"
        report["error_type"] = type(error).__name__
        report["error"] = str(error)
        print(f"REVISION_VALIDATION_FAILED: {type(error).__name__}: {error}", file=sys.stderr)
    finally:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if report["status"] != "passed":
        print(f"Audit JSON: {report_path}", file=sys.stderr)
        return 1
    print("REVISION_VALIDATION_PASSED: 9 target scales 0-2; 9 other scales/PNGs unchanged; "
          "6 RDS files unchanged; 26 PDFs / 31 pages rendered at 100 dpi")
    print(f"Audit JSON: {report_path}")
    print(f"Contact sheet: {report['contact_sheet']['status']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
