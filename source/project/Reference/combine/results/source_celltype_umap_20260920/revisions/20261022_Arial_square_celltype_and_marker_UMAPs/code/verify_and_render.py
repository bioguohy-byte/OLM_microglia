#!/usr/bin/env python3
"""Verify PDF dimensions/Arial/inputs, render every page for visual QA."""
import csv, hashlib, json, math, pathlib, re, subprocess, sys
import xml.etree.ElementTree as ET
from PIL import Image, ImageDraw, ImageOps
task = pathlib.Path(sys.argv[1]).resolve()
audit = task / 'audit'
render = audit / 'pdf_renders'
render.mkdir(exist_ok=True)
def rows(name):
    with (audit / name).open(newline='') as fh:
        return list(csv.DictReader(fh))
def sha(path): return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()
panel = rows('physical_panel_dimensions.csv')
assert len(panel) == 24
assert {r['font'] for r in panel} == {'Arial'}
assert all(math.isclose(float(r[k]),114.3,abs_tol=1e-7) for r in panel for k in ('panel_width_mm','panel_height_mm'))
bars = rows('selected_marker_bar_provenance.csv')
scales = rows('marker_scale_parameters.csv')
assert len(bars)==len(scales)==18
for row in bars:
    assert sha(task / 'inputs/selected_marker_pdf' / row['zip_filename']) == row['sha256'] == row['historical_sha256']
    scale = next(r for r in scales if r['study']==row['study'] and r['marker']==row['marker'])
    for key in ('lower_limit','upper_limit'):
        assert math.isclose(float(scale[key]),float(row[key]),rel_tol=1e-13)
    assert scale['colour_low']==row['colour_low'] and scale['colour_high']==row['colour_high']
    assert scale['visible_breaks']==row['tick_values']
    assert scale['visible_tick_labels']==row['tick_labels']
for row in rows('input_copy_sha256.csv'):
    assert sha(task / row['destination']) == row['sha256'] == sha(row['source'])
coordinate = rows('coordinate_and_annotation_preservation.csv')
assert len(coordinate)==6
for row in coordinate:
    assert row['umap_recomputed']=='FALSE' and row['annotation_modified']=='FALSE'
    assert row['celltype_and_marker_coordinates_identical_by_barcode']=='TRUE'
    assert math.isclose(float(row['view_x_max'])-float(row['view_x_min']),
                        float(row['view_y_max'])-float(row['view_y_min']),rel_tol=1e-12)
manifest = rows('output_manifest.csv')
assert len(manifest)==27
assert len(list((task / 'output/pdf').glob('*.pdf')))==27
assert len(list((task / 'output/png').glob('*.png')))==27
pdf_audit=[]
for row in manifest:
    pdf = task / row['pdf']
    info = subprocess.check_output(['pdfinfo',str(pdf)],text=True)
    assert int(re.search(r'^Pages:\s+(\d+)',info,re.M)[1])==1
    width,height = map(float,re.search(r'^Page size:\s+([\d.]+) x ([\d.]+)',info,re.M).groups())
    # R's Cairo 1.16 device truncates the outer PDF canvas to whole points.
    # The fixed grid panel remains exactly 324 points (114.3 mm) per side.
    assert math.isclose(width,math.floor(float(row['width_in'])*72),abs_tol=.03)
    assert math.isclose(height,math.floor(float(row['height_in'])*72),abs_tol=.03)
    fonts = subprocess.check_output(['pdffonts',str(pdf)],text=True)
    lines=fonts.splitlines()[2:]
    assert lines and all('Arial' in line for line in lines),fonts
    assert all(re.search(r'\byes\s+yes\b',line) for line in lines),fonts
    (audit / (pdf.name+'.fonts.txt')).write_text(fonts)
    text = subprocess.check_output(['pdftotext','-layout',str(pdf),'-'],text=True)
    assert 'UMAP_1' in text and 'UMAP_2' in text
    (audit / (pdf.name+'.text.txt')).write_text(text)
    if row['name'].endswith('_FeaturePlot'):
        study,gene,_=row['name'].split('_',2)
        chosen=next(r for r in bars if r['study']==study and r['marker']==gene)
        bbox=subprocess.check_output(['pdftotext','-bbox',str(pdf),'-'],text=True)
        words=list(ET.fromstring(bbox).iter('{http://www.w3.org/1999/xhtml}word'))
        rna=next(w for w in words if w.text=='RNA')
        rx=float(rna.attrib['xMin']); ry=float(rna.attrib['yMax'])
        ticks=[w for w in words if re.fullmatch(r'\d+(?:\.\d+)?',w.text or '')
               and float(w.attrib['xMin'])>=rx-2 and ry<float(w.attrib['yMin'])<ry+120]
        ticks.sort(key=lambda w:float(w.attrib['yMin']),reverse=True)
        assert ';'.join(w.text for w in ticks)==chosen['tick_labels'],(pdf.name,[w.text for w in ticks])
    target = render / row['name']
    subprocess.run(['pdftoppm','-f','1','-l','1','-scale-to','2400' if int(row['n_panels'])>1 else '1500',
                    '-singlefile','-png',str(pdf),str(target)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    pdf_audit.append(dict(file=row['pdf'],sha256=sha(pdf),pages=1,width_pt=width,height_pt=height,
                          all_fonts_Arial=True,fonts_embedded=True,pdf_rendered=True))
with (audit/'pdf_validation.csv').open('w',newline='') as fh:
    writer=csv.DictWriter(fh,fieldnames=list(pdf_audit[0]));writer.writeheader();writer.writerows(pdf_audit)

# A compact index of all 18 selected originals complements the exact hash audit.
tiles=[]
for study in ('Chia','Fullard','Grubman','Jakel','Su','Yao'):
    for gene in ('P2RY12','CX3CR1','CSF1R'):
        row=next(r for r in bars if r['study']==study and r['marker']==gene)
        dest=render/('selected_'+study+'_'+gene)
        subprocess.run(['pdftoppm','-f','1','-l','1','-scale-to','700','-singlefile','-png',
                        str(task/'inputs/selected_marker_pdf'/row['zip_filename']),str(dest)],check=True,
                        stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        image=Image.open(str(dest)+'.png').convert('RGB')
        tile=Image.new('RGB',(740,750),'white')
        thumb=ImageOps.contain(image,(720,710))
        tile.paste(thumb,((740-thumb.width)//2,25))
        ImageDraw.Draw(tile).text((12,8),study+' / '+gene+' / 0-'+row['upper_limit'],fill='black')
        tiles.append(tile)
sheet=Image.new('RGB',(3*740,6*750),'#eeeeee')
for i,tile in enumerate(tiles): sheet.paste(tile,((i%3)*740,(i//3)*750))
sheet.save(audit/'selected_marker_originals_contact_sheet.png')
report=dict(status='PASS',selected_source_PDF_exact_hash_matches=18,display_scales_preserved=18,
            coordinates_and_annotation_preserved_studies=6,identical_square_data_panels=24,
            panel_side_mm=114.3,output_PDFs_verified_and_rendered=27,all_fonts='embedded Arial',
            celltype_grid='3 columns x 2 rows',marker_grids='2 PDFs; 3 columns x 3 rows each')
(audit/'verification_summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
