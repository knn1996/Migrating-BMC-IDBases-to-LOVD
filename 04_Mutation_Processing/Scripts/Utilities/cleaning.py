import zipfile, shutil, os

base = r"C:\Users\BornLoser\Desktop\pipeline_doc_v6_extracted"
out  = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\06_Writing\pipeline_documentation_v7.docx"
with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zf:
    for root, dirs, files in os.walk(base):
        for file in files:
            fp = os.path.join(root, file)
            zf.write(fp, os.path.relpath(fp, base))
print("Packed to", out)