# -*- mode: python ; coding: utf-8 -*-

from PyInstaller.utils.hooks import collect_data_files, collect_submodules


hiddenimports = collect_submodules("Bio") + [
    "certifi", "openpyxl", "PIL.Image", "PIL.ImageTk", "screeninfo",
    "tkinter", "tkinter.filedialog", "tkinter.messagebox", "tkinter.ttk",
    "tqdm.auto", "yaml",
]

datas = collect_data_files("Bio") + [
    ("pynnotate/data/pynnotate.png", "data"),
    ("pynnotate/data/synonyms_database.json", "data"),
    ("pynnotate/data/synonyms_manifest.json", "data"),
    ("pynnotate/examples", "examples"),
]

a = Analysis(
    ["pynnotate/pynnotate_graphic.py"],
    pathex=["."], binaries=[], datas=datas, hiddenimports=hiddenimports,
    hookspath=[], hooksconfig={}, runtime_hooks=[],
    excludes=["Bio.PDB.mmtf", "mmtf", "Bio.Alphabet"], noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz, a.scripts, [], exclude_binaries=True, name="Pynnotate",
    debug=False, bootloader_ignore_signals=False, strip=False, upx=False,
    console=False, icon="pynnotate/data/pynnotate.png",
)

coll = COLLECT(
    exe, a.binaries, a.datas, strip=False, upx=False, name="Pynnotate",
)
