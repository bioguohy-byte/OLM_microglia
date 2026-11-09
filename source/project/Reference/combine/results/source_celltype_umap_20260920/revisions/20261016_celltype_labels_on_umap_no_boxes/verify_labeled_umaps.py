#!/usr/bin/python3
"""Validate the saved labelled UMAP figures; never alter inputs or plot data."""

import argparse
import csv
import datetime
import json
import math
import re
import shutil
import subprocess
import sys
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path


MAIN_STUDIES = {"Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao"}
BOOKS = {
    "main": "Six_study_celltypes_labeled_on_UMAPs.pdf",
    "supplementary": "Supplementary_existing_annotation_columns_on_UMAPs.pdf",
}


class Validator:
    def __init__(self, root):
        self.root = root
        self.report = {
            "task_dir": str(root),
            "checked_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "checks": [],
            "errors": [],
            "pdfs": [],
        }

    def check(self, condition, description, **details):
        ok = bool(condition)
        self.report["checks"].append(dict(check=description, passed=ok, **details))
        if not ok:
            suffix = ": " + json.dumps(details, ensure_ascii=False) if details else ""
            self.report["errors"].append(description + suffix)
        return ok

    def table(self, name, columns, delimiter=","):
        path = self.root / name
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            missing = sorted(set(columns) - set(reader.fieldnames or []))
            if missing:
                raise ValueError("{} is missing columns: {}".format(name, ", ".join(missing)))
            return list(reader)

    def path(self, relative):
        path = (self.root / relative).resolve()
        if self.root != path and self.root not in path.parents:
            raise ValueError("Path escapes the task directory: " + relative)
        return path

    def command(self, args):
        result = subprocess.run(args, capture_output=True, text=True, timeout=180)
        if result.returncode:
            raise RuntimeError("{} failed: {}".format(args[0], result.stderr.strip()))
        return result.stdout

    def pdf(self, path, page_count, width_in, height_in, render_name):
        relative = str(path.relative_to(self.root))
        record = {"file": relative, "expected_pages": page_count}
        self.report["pdfs"].append(record)
        if not self.check(path.is_file(), "PDF exists: " + relative):
            return None, []
        try:
            info = self.command(["pdfinfo", str(path)])
            match = re.search(r"(?m)^Pages:\s*(\d+)", info)
            if not match:
                raise ValueError("pdfinfo did not report a page count")
            actual_pages = int(match.group(1))
            record["actual_pages"] = actual_pages
            self.check(actual_pages == page_count, "PDF page count: " + relative,
                       expected=page_count, actual=actual_pages)
            box_info = self.command(["pdfinfo", "-f", "1", "-l", str(actual_pages),
                                     "-box", str(path)])
            sizes = [(float(w), float(h)) for w, h in re.findall(
                r"(?m)^Page(?:\s+\d+)?\s+size:\s*([\d.]+)\s+x\s+([\d.]+)\s+pts", box_info)]
            expected = [width_in * 72, height_in * 72]
            record["page_sizes_pt"] = sizes
            record["expected_size_pt"] = expected
            self.check(len(sizes) == actual_pages, "PDF reports every page size: " + relative,
                       expected=actual_pages, actual=len(sizes))
            self.check(bool(sizes) and all(abs(w - expected[0]) < 1 and
                                          abs(h - expected[1]) < 1 for w, h in sizes),
                       "PDF dimensions within Cairo whole-point tolerance (<1 pt): " + relative,
                       expected=expected, actual=sizes)
            text = self.command(["pdftotext", "-enc", "UTF-8", "-layout", str(path), "-"])
            # Spatial extraction interleaves lines from neighbouring multi-line
            # labels. Raw drawing order preserves the text within each label.
            text += "\n" + self.command(["pdftotext", "-enc", "UTF-8", "-raw", str(path), "-"])
            destination = self.root / "audit" / "pdf_render" / render_name
            destination.mkdir(parents=True, exist_ok=True)
            # Remove only this verifier's previous page images before rendering again.
            for old in destination.glob("page-*.png"):
                old.unlink()
            self.command(["pdftoppm", "-r", "100", "-png", str(path), str(destination / "page")])
            pages = sorted(destination.glob("page-*.png"), key=lambda p: int(p.stem.rsplit("-", 1)[1]))
            self.check(len(pages) == actual_pages, "Rendered all PDF pages at 100 dpi: " + relative,
                       expected=actual_pages, actual=len(pages))
            record["rendered_pages"] = [str(p.relative_to(self.root)) for p in pages]
            return text, pages
        except Exception as exc:
            self.check(False, "PDF inspection/render failed: " + relative, reason=str(exc))
            return None, []

    def run(self):
        from PIL import Image

        for command in ("pdfinfo", "pdftotext", "pdftoppm"):
            if not shutil.which(command):
                raise RuntimeError("Required existing tool is unavailable: " + command)
        manifest = self.table("plot_manifest.tsv", ["study", "bundle_file", "annotation_column",
                              "display_group", "plot_id", "title_note"], "\t")
        dims = self.table("audit/plot_dimensions.csv", ["plot_id", "pdf", "png", "width_in",
                          "height_in", "png_dpi", "n_cells", "n_labels"])
        counts = self.table("tables/celltype_counts_and_colours.csv", ["plot_id", "study",
                            "display_group", "annotation_column", "celltype", "n_cells", "colour"])
        anchors = self.table("tables/label_anchors_and_display.csv", ["celltype", "n_cells", "x", "y",
                             "anchor_cell", "label_text", "nudge_x", "nudge_y", "colour", "plot_id",
                             "annotation_column"])
        preservation = self.table("audit/plot_preservation.tsv", ["plot_id", "study", "annotation_column",
                                   "n_cells", "n_labels", "coordinate_ratio", "pointsize",
                                   "coordinates_identical_to_bundle", "labels_identical_to_bundle",
                                   "umap_recomputed", "labels_reassigned", "legend_shown"], "\t")
        ids = [row["plot_id"] for row in manifest]
        self.check(bool(ids) and len(ids) == len(set(ids)), "Manifest has unique nonempty plot IDs")
        self.check(all(row["display_group"] in BOOKS for row in manifest),
                   "Manifest display_group values are main or supplementary")
        main = [row for row in manifest if row["display_group"] == "main"]
        supplementary = [row for row in manifest if row["display_group"] == "supplementary"]
        self.check(len(main) == 6 and {row["study"] for row in main} == MAIN_STUDIES,
                   "Exactly six main studies", expected=sorted(MAIN_STUDIES),
                   actual=[row["study"] for row in main])
        self.check(bool(supplementary), "Supplementary annotation panels are listed")
        expected_ids = set(ids)
        for label, rows in (("dimensions", dims), ("counts", counts),
                            ("anchors", anchors), ("preservation", preservation)):
            self.check({row["plot_id"] for row in rows} == expected_ids,
                       "Complete plot coverage in " + label,
                       missing=sorted(expected_ids - {row["plot_id"] for row in rows}),
                       unexpected=sorted({row["plot_id"] for row in rows} - expected_ids))
        for label, rows in (("dimensions", dims), ("preservation", preservation)):
            self.check(all(n == 1 for n in Counter(row["plot_id"] for row in rows).values()),
                       "One row per plot in " + label)

        grouped_counts, grouped_anchors = defaultdict(list), defaultdict(list)
        for row in counts:
            grouped_counts[row["plot_id"]].append(row)
        for row in anchors:
            grouped_anchors[row["plot_id"]].append(row)
        by_dims = {row["plot_id"]: row for row in dims}
        by_preservation = {row["plot_id"]: row for row in preservation}
        expected_pdf_paths, main_renders = set(), {}
        for row in manifest:
            plot_id = row["plot_id"]
            if plot_id not in by_dims or plot_id not in by_preservation:
                continue
            dim, kept = by_dims[plot_id], by_preservation[plot_id]
            type_rows, label_rows = grouped_counts[plot_id], grouped_anchors[plot_id]
            n_cells, n_labels = integer(dim["n_cells"]), integer(dim["n_labels"])
            width, height, dpi = float(dim["width_in"]), float(dim["height_in"]), float(dim["png_dpi"])
            self.check(n_cells > 0 and n_labels > 0, plot_id + ": positive cell and label counts")
            self.check(width > 0 and height > 0 and dpi == 240, plot_id + ": dimensions and 240 dpi",
                       width_in=width, height_in=height, dpi=dpi)
            self.check(self.path("plot_data/" + row["bundle_file"]).is_file(),
                       plot_id + ": source plot bundle exists")
            self.check(len(type_rows) == n_labels and
                       len({r["celltype"] for r in type_rows}) == n_labels,
                       plot_id + ": exactly one count row per cell type", expected=n_labels,
                       actual=len(type_rows))
            self.check(all(integer(r["n_cells"]) > 0 for r in type_rows) and
                       sum(integer(r["n_cells"]) for r in type_rows) == n_cells,
                       plot_id + ": cell-type counts sum to the plotted cell count", expected=n_cells,
                       actual=sum(integer(r["n_cells"]) for r in type_rows))
            self.check(len(label_rows) == n_labels and
                       Counter(r["celltype"] for r in label_rows) == Counter(r["celltype"] for r in type_rows),
                       plot_id + ": one anchor for every cell type without omissions", expected=n_labels,
                       actual=len(label_rows))
            self.check(len({r["anchor_cell"] for r in label_rows}) == len(label_rows) and
                       all(r["anchor_cell"] for r in label_rows), plot_id + ": anchor cells are unique and named")
            self.check(all(math.isfinite(float(r[key])) for r in label_rows
                           for key in ("x", "y", "nudge_x", "nudge_y")),
                       plot_id + ": label anchor coordinates and nudges are finite")
            self.check(all(r["annotation_column"] == row["annotation_column"]
                           for r in type_rows + label_rows + [kept]) and
                       all(r["study"] == row["study"] and r["display_group"] == row["display_group"]
                           for r in type_rows) and kept["study"] == row["study"],
                       plot_id + ": study and annotation column agree with the manifest")
            type_map = {r["celltype"]: r for r in type_rows}
            self.check(all(r["celltype"] in type_map and
                           integer(r["n_cells"]) == integer(type_map[r["celltype"]]["n_cells"]) and
                           r["colour"] == type_map[r["celltype"]]["colour"] and
                           compact(r["celltype"]) in compact(r["label_text"]) and
                           count_in_text(r["label_text"], integer(r["n_cells"])) for r in label_rows),
                       plot_id + ": displayed labels preserve type names, counts, and colours")
            for flag, expected in (("coordinates_identical_to_bundle", True),
                                   ("labels_identical_to_bundle", True), ("umap_recomputed", False),
                                   ("labels_reassigned", False), ("legend_shown", False)):
                self.check(boolean(kept[flag]) == expected, plot_id + ": " + flag,
                           expected=expected, actual=kept[flag])
            self.check(math.isclose(float(kept["coordinate_ratio"]), 1, abs_tol=1e-12) and
                       math.isclose(float(kept["pointsize"]), 2.2, abs_tol=1e-12),
                       plot_id + ": coordinate ratio 1 and point size 2.2")
            self.check(integer(kept["n_cells"]) == n_cells and integer(kept["n_labels"]) == n_labels,
                       plot_id + ": preservation audit agrees with figure counts")

            pdf_path, png_path = self.path(dim["pdf"]), self.path(dim["png"])
            expected_pdf_paths.add(pdf_path)
            text, pages = self.pdf(pdf_path, 1, width, height, "single_" + plot_id)
            if row["display_group"] == "main" and pages:
                main_renders[plot_id] = pages[0]
            if text is not None:
                absent_names = [r["celltype"] for r in type_rows if compact(r["celltype"]) not in compact(text)]
                absent_counts = [integer(r["n_cells"]) for r in type_rows
                                 if not count_in_text(text, integer(r["n_cells"]))]
                self.check(not absent_names, plot_id + ": all cell-type names occur in PDF text", missing=absent_names)
                self.check(not absent_counts, plot_id + ": all parenthesized cell counts occur in PDF text",
                           missing=absent_counts)
            if self.check(png_path.is_file(), "PNG exists: " + dim["png"]):
                with Image.open(png_path) as image:
                    expected_pixels = (round(width * 240), round(height * 240))
                    self.check(image.size == expected_pixels, plot_id + ": PNG pixel dimensions at 240 dpi",
                               expected=expected_pixels, actual=image.size)
                    saved_dpi = image.info.get("dpi")
                    self.check(saved_dpi is not None and all(abs(x - 240) < 0.1 for x in saved_dpi),
                               plot_id + ": PNG resolution metadata is 240 dpi", actual=saved_dpi)
                    image.verify()

        for group, filename in BOOKS.items():
            book_rows = [r for r in manifest if r["display_group"] == group]
            if not book_rows:
                continue
            path = self.root / filename
            expected_pdf_paths.add(path)
            book_dims = [by_dims[r["plot_id"]] for r in book_rows if r["plot_id"] in by_dims]
            self.check(len(book_dims) == len(book_rows) and len({(d["width_in"], d["height_in"])
                       for d in book_dims}) == 1, group + ": book panels have consistent page dimensions")
            if book_dims:
                self.pdf(path, len(book_rows), float(book_dims[0]["width_in"]),
                         float(book_dims[0]["height_in"]), "book_" + path.stem)
        actual_pdf_paths = {p.resolve() for p in self.root.rglob("*.pdf") if p.is_file()}
        self.check(actual_pdf_paths == expected_pdf_paths, "All PDFs are accounted for and rendered",
                   unexpected=sorted(str(p.relative_to(self.root)) for p in actual_pdf_paths - expected_pdf_paths),
                   missing=sorted(str(p.relative_to(self.root)) for p in expected_pdf_paths - actual_pdf_paths))
        self.report["manifest_summary"] = {"panels": len(manifest), "main_panels": len(main),
                                           "supplementary_panels": len(supplementary),
                                           "expected_pdf_files": len(expected_pdf_paths),
                                           "expected_pdf_pages": 2 * len(manifest)}
        self.contact_sheet(main, main_renders)

    def contact_sheet(self, main, renders):
        from PIL import Image, ImageDraw, ImageFont

        if not self.check(len(main) == 6 and all(r["plot_id"] in renders for r in main),
                          "Six main PDF renders are available for the contact sheet"):
            return
        tile_w, tile_h, gap = 800, 680, 16
        sheet = Image.new("RGB", (2 * tile_w + 3 * gap, 3 * tile_h + 4 * gap), "#e8e8e8")
        draw = ImageDraw.Draw(sheet)
        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 28)
        except OSError:
            font = ImageFont.load_default()
        resample = getattr(Image, "Resampling", Image).LANCZOS
        for index, row in enumerate(main):
            x, y = gap + (index % 2) * (tile_w + gap), gap + (index // 2) * (tile_h + gap)
            draw.rectangle((x, y, x + tile_w, y + tile_h), fill="white")
            draw.text((x + 18, y + 10), row["study"], fill="black", font=font)
            with Image.open(renders[row["plot_id"]]) as page:
                thumbnail = page.convert("RGB")
                thumbnail.thumbnail((tile_w - 24, tile_h - 55), resample)
                sheet.paste(thumbnail, (x + (tile_w - thumbnail.width) // 2,
                                       y + 48 + (tile_h - 55 - thumbnail.height) // 2))
        destination = self.root / "audit" / "contact_sheet.png"
        sheet.save(destination)
        self.report["contact_sheet"] = str(destination.relative_to(self.root))

    def save(self):
        self.report["passed"] = not self.report["errors"]
        destination = self.root / "audit" / "figure_validation.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(self.report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        return destination


def compact(text):
    return re.sub(r"\s+", "", unicodedata.normalize("NFKC", text))


def integer(value):
    number = float(value)
    if not math.isfinite(number) or not number.is_integer():
        raise ValueError("Expected an integer, got: " + str(value))
    return int(number)


def boolean(value):
    normalized = value.strip().lower()
    if normalized not in ("true", "false", "1", "0"):
        raise ValueError("Expected an audit boolean, got: " + str(value))
    return normalized in ("true", "1")


def count_in_text(text, n_cells):
    normalized = compact(text).replace(",", "")
    return "({})".format(n_cells) in normalized


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_dir", nargs="?", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    root = args.task_dir.resolve(strict=True)
    validator = Validator(root)
    try:
        validator.run()
    except Exception as exc:
        validator.check(False, "Validation could not finish", reason="{}: {}".format(type(exc).__name__, exc))
    report_path = validator.save()
    if validator.report["passed"]:
        summary = validator.report["manifest_summary"]
        print("PASS: {panels} panels, {expected_pdf_files} PDFs, {expected_pdf_pages} pages; "
              "all dimensions, labels, counts and preservation checks passed.".format(**summary))
    else:
        print("FAIL: {} validation problem(s):".format(len(validator.report["errors"])), file=sys.stderr)
        for error in validator.report["errors"]:
            print("- " + error, file=sys.stderr)
    print("Validation report: " + str(report_path))
    return 0 if validator.report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
