# -*- mode: python ; coding: utf-8 -*-
# PyInstaller 打包配置：基因型分析 GUI -> Windows exe

block_cipher = None

# 打包 parasail 的 DLL，否则 exe 内 parasail 无法加载，会走回退比对导致结果与直接运行 Python 不一致
from PyInstaller.utils.hooks import collect_all
_parasail_datas, parasail_binaries, _parasail_hidden = collect_all('parasail')

# 打包时需包含的数据文件（不依赖 known_polymorphisms.txt，程序内无该文件时使用空列表）
added_files = []

# 排除无关大库，减小 exe 体积（Bio/parasail/numpy 保留）
excludes_list = [
    'pandas', 'scipy', 'matplotlib', 'PIL', 'cv2', 'skimage',
    'torch', 'tensorflow', 'keras', 'transformers', 'dask', 'bokeh',
    'IPython', 'jupyter', 'notebook', 'pytest', 'sphinx', 'docutils',
    'PyQt5', 'PyQt6', 'PySide2', 'PySide6', 'wx', 'zmq',
    'xlrd', 'openpyxl', 'xlsxwriter', 'pyarrow', 'tables',
    'Crypto', 'cryptography', 'nacl', 'bcrypt', 'paramiko',
    'h5py', 'netCDF4', 'pyproj', 'geopandas',
    'numba', 'llvmlite', 'cython', 'setuptools',
    'pip', 'wheel', 'pkg_resources', 'importlib_metadata',
]

a = Analysis(
    ['genotype_gui.py'],
    pathex=[],
    binaries=parasail_binaries,
    datas=added_files,
    hiddenimports=[
        'Bio.SeqIO',
        'Bio.Align',
        'Bio.pairwise2',
        'parasail',
        'pysanger',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=excludes_list,
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='GeneType',
    icon=None,  # 不嵌入图标，避免杀毒软件在写入资源时报错
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,   # 不显示控制台窗口（GUI 程序）
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
