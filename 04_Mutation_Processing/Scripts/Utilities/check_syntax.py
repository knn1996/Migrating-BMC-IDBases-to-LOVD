import py_compile, os

base = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Scripts"
scripts = []
for root, dirs, files in os.walk(base):
    for f in files:
        if f.endswith(".py") and f != "check_syntax.py":
            scripts.append(os.path.join(root, f))

scripts.sort()
ok, errors = [], []
for s in scripts:
    rel = s.replace(base + "\\", "")
    try:
        py_compile.compile(s, doraise=True)
        ok.append(rel)
        print(f"OK      {rel}")
    except py_compile.PyCompileError as e:
        errors.append((rel, str(e)))
        print(f"ERROR   {rel}\n        {e}")

print(f"\n--- {len(ok)} OK, {len(errors)} ERROR ---")
