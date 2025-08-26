import sys
from pathlib import Path

# Ensure the ChEMBL package is importable during tests
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "ChEMBL"))
