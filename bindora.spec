# -*- mode: python ; coding: utf-8 -*-
import sys
import os
from pathlib import Path
from PyInstaller.utils.hooks import collect_all, collect_data_files, collect_submodules

block_cipher = None

# Collect complex scientific packages
rdkit_datas, rdkit_binaries, rdkit_hidden = collect_all('rdkit')
meeko_datas, meeko_binaries, meeko_hidden = collect_all('meeko')
gemmi_datas, gemmi_binaries, gemmi_hidden = collect_all('gemmi')
waitress_datas, waitress_binaries, waitress_hidden = collect_all('waitress')
webview_datas, webview_binaries, webview_hidden = collect_all('webview')

# Collect OpenMM package and runtime DLL/plugin libraries
try:
    openmm_datas, openmm_binaries, openmm_hidden = collect_all('openmm')
except Exception:
    openmm_datas, openmm_binaries, openmm_hidden = [], [], []

openmm_extra_datas = []
try:
    import openmm
    openmm_root = os.path.dirname(openmm.__file__)
    openmm_libs_path = os.path.abspath(os.path.join(openmm_root, "..", "OpenMM.libs"))
    if os.path.exists(openmm_libs_path):
        openmm_extra_datas.append((openmm_libs_path, "OpenMM.libs"))
except Exception:
    pass

# Application assets
datas = [
    ('frontend', 'frontend'),
    ('bin/vina.exe', 'bin'),
    ('backend/utils/fpscores.pkl.gz', 'backend/utils'),
    ('backend/services/boiled_egg_coords.json', 'backend/services'),
    ('data/benchmarks', 'data/benchmarks'),
]

datas += rdkit_datas + meeko_datas + gemmi_datas + openmm_datas + openmm_extra_datas + waitress_datas + webview_datas
binaries = rdkit_binaries + meeko_binaries + gemmi_binaries + openmm_binaries + waitress_binaries + webview_binaries

hiddenimports = [
    'waitress',
    'webview',
    'clr',
    'pythonnet',
    'openmm',
    'openmm.app',
    'openmm.unit',
    'scipy',
    'scipy.spatial',
    'scipy.spatial.distance',
    'scipy.special',
    'sklearn',
    'sklearn.utils._typedefs',
    'sqlalchemy',
    'sqlalchemy.sql.default_comparator',
    'sqlalchemy.dialects.sqlite',
    'flask',
    'flask_cors',
] + rdkit_hidden + meeko_hidden + gemmi_hidden + openmm_hidden + waitress_hidden + webview_hidden

a = Analysis(
    ['desktop_launcher.py'],
    pathex=['.'],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['pytest', 'IPython', 'jupyter', 'sphinx'],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='bindora_launcher',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='frontend/assets/branding/favicon.ico',
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='bindora_launcher',
)
