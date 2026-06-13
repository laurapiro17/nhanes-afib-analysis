"""Download the NHANES 2017-2018 public data files used in this analysis.

Run once; files are cached in ../data and skipped if already present. These are
public-domain U.S. government data (CDC / NCHS), redistributed here only as a
download script rather than bundling the raw files.
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles"
FILES = {
    "DEMO_J.xpt": "Demographics (age, sex, race/ethnicity, survey weights & design)",
    "BPX_J.xpt": "Blood pressure examination (measured SBP/DBP)",
    "BMX_J.xpt": "Body measures (BMI)",
    "BPQ_J.xpt": "Blood pressure questionnaire (antihypertensive medication use)",
}
DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def main() -> int:
    DATA_DIR.mkdir(exist_ok=True)
    for name, desc in FILES.items():
        dest = DATA_DIR / name
        if dest.exists() and dest.stat().st_size > 100_000:
            print(f"[skip] {name} already present ({dest.stat().st_size:,} bytes)")
            continue
        url = f"{BASE}/{name}"
        print(f"[get ] {name} <- {url}  ({desc})")
        urllib.request.urlretrieve(url, dest)
        size = dest.stat().st_size
        if size < 100_000:
            print(f"[err ] {name} is only {size} bytes - likely not the data file", file=sys.stderr)
            return 1
        print(f"       {size:,} bytes")
    print("All files ready in", DATA_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
