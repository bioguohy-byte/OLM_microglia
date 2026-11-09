#!/usr/bin/env python3
import csv, hashlib, json, pathlib, sys, zipfile
task=pathlib.Path(sys.argv[1]).resolve()
assert json.loads((task/'audit/verification_summary.json').read_text())['status']=='PASS'
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
files=[]
for folder in ('inputs','code','output','audit','runtime/fonts'):
    for p in sorted((task/folder).rglob('*')):
        rel=str(p.relative_to(task))
        if not p.is_file(): continue
        if rel.startswith(('audit/pdf_renders/','audit/selected_bar_render/','audit/selected_bar_bbox/')): continue
        if p.suffix=='.png' and p.name!='selected_marker_originals_contact_sheet.png' and folder=='audit': continue
        if p.name in ('deliverable_file_sha256.csv','package_sha256.txt'): continue
        files.append(p)
files.append(task/'README.txt')
inventory=task/'audit/deliverable_file_sha256.csv'
with inventory.open('w',newline='') as fh:
    writer=csv.DictWriter(fh,fieldnames=['path','bytes','sha256'])
    writer.writeheader()
    writer.writerows(dict(path=str(p.relative_to(task)),bytes=p.stat().st_size,sha256=sha(p)) for p in files)
files.append(inventory)
archive=task/'Six_study_UMAPs_Arial_square_20261022.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files: z.write(p,str(p.relative_to(task)))
digest=sha(archive)
(task/'audit/package_sha256.txt').write_text(digest+'  '+archive.name+'\n')
print(json.dumps(dict(archive=str(archive),bytes=archive.stat().st_size,sha256=digest,files=len(files)),indent=2))
