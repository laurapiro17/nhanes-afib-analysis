import sys
from pathlib import Path

# The analysis is a plain script under src/ (no installed package); put it on the
# path so tests can import its functions directly.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
