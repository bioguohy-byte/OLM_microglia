#!/usr/bin/env python3
"""Copy selected thin plotting caches and audit every selected PDF by SHA256."""
import csv, hashlib, pathlib, shutil, subprocess, sys, zipfile

task = pathlib.Path(sys.argv[1]).resolve()
combine = pathlib.Path(sys.argv[2]).resolve()
base = combine / 'results/source_celltype_umap_20260920/revisions/20261016_celltype_labels_on_umap_no_boxes'
yao = base.parent / '20261016_Yao_author_celltypes'
marker = combine / 'results/six_study_microglia_marker_FeaturePlot_20261014'
for name in ('inputs/celltype', 'inputs/marker', 'inputs/selected_marker_pdf', 'inputs/source_tables', 'audit'):
    (task / name).mkdir(parents=True, exist_ok=True)
copies = []
def copy(src, dst):
    shutil.copy2(src, dst)
    copies.append(dict(source=str(src), destination=str(dst.relative_to(task)), sha256=sha(src)))
def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()
studies = ['Chia', 'Fullard', 'Grubman', 'Jakel', 'Su', 'Yao']
for study in studies:
    bundle = ('Jakel_featureplot_aligned_labels_umap.rds' if study == 'Jakel' else study + '_labels_umap.rds')
    src = yao / 'plot_data/Yao_current_UMAP_author_labels.rds' if study == 'Yao' else base / 'plot_data' / bundle
    copy(src, task / 'inputs/celltype' / (study + '.rds'))
    copy(marker / 'plot_data' / (study + '_stored_umap_three_RNA_markers.rds'), task / 'inputs/marker' / (study + '.rds'))
for src, dest in [(base / 'tables/celltype_counts_and_colours.csv', 'base_celltype_colours.csv'),
                  (base / 'tables/label_anchors_and_display.csv', 'base_celltype_anchors.csv'),
                  (yao / 'audit/plot_counts_and_colours.csv', 'yao_celltype_colours.csv'),
                  (yao / 'audit/plot_label_anchors.csv', 'yao_celltype_anchors.csv'),
                  (yao / 'code/label_position_overrides.csv', 'yao_label_overrides.csv'),
                  (marker / 'audit/marker_expression_ranges.tsv', 'baseline_marker_expression_ranges.tsv')]:
    copy(src, task / 'inputs/source_tables' / dest)

upper = dict(P2RY12=4.2219674636854, CX3CR1=4.45026076501754, CSF1R=4.35729671715264)
audit = []
with zipfile.ZipFile(task / 'inputs/selected_marker_figures.zip') as archive:
    names = [n for n in archive.namelist() if n.lower().endswith('.pdf') and not n.startswith('__MACOSX/')]
    assert len(names) == 18, names
    for name in names:
        filename = pathlib.PurePosixPath(name).name
        study, gene, _ = filename.split('_', 2)
        assert study in studies and gene in upper
        destination = task / 'inputs/selected_marker_pdf' / filename
        destination.write_bytes(archive.read(name))
        if study in ('Grubman', 'Jakel'):
            history = marker / 'revisions/20261016_Grubman_Jakel_range0_1_vs_1p5/range0_1'
            cap = 1.0
        elif study == 'Fullard':
            history = marker / 'revisions/20261016_Grubman_Jakel_Fullard_range0_2'
            cap = 2.0
        else:
            history = marker
            cap = upper[gene]
        history = history / 'pdf_single_marker' / (study + '_' + gene + '_FeaturePlot.pdf')
        digest = sha(destination)
        assert digest == sha(history), (filename, history)
        pdf_text = subprocess.check_output(['pdftotext', '-layout', str(destination), '-'], text=True)
        (task / 'audit' / (filename + '.txt')).write_text(pdf_text)
        tick_values = '0;0.5;1' if cap == 1 else ('0;1;2' if cap == 2 else '0;2;4')
        tick_labels = '0.0;0.5;1.0' if cap == 1 else tick_values
        audit.append(dict(study=study, marker=gene, zip_filename=filename, sha256=digest,
                          historical_file=str(history), historical_sha256=sha(history), exact_match=True,
                          lower_limit=0, upper_limit=cap, tick_values=tick_values,tick_labels=tick_labels,
                          colour_low='#D3D3D3', colour_high='#E31A1C'))
assert len({(r['study'],r['marker']) for r in audit}) == 18
for name, rows in [('selected_marker_bar_provenance.csv', audit), ('input_copy_sha256.csv', copies)]:
    with (task / 'audit' / name).open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
print('PREPARED: 12 thin bundles; 18/18 selected marker PDF hashes match historical exports')
