"""Generate a small, synthetic two-pool mzXML example without experimental data."""

import argparse
import base64
import struct
from pathlib import Path

import numpy as np
import pandas as pd
from pyteomics import mass


def generate_example(out, pools=("P1", "P2")):
    """Create Gaussian MS1 traces: standards and T1 are stronger than T10."""
    directory = Path(out)
    raw = directory / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    rows = []
    for pool, formula, smiles, carbons in zip(
        pools, ["C6H6", "C7H8"], ["c1ccccc1", "Cc1ccccc1"], [6, 7]
    ):
        neutral = mass.calculate_mass(formula=formula)
        target = neutral + 1.007285
        rows.append(
            {
                "SGC ID for Component": f"CMP_{pool}",
                "SGC ID for Pool": pool,
                "Comp Acc Mass": neutral,
                "Pooled Well": pool,
                "SMILES": smiles,
                "formula": formula,
            }
        )
        for group, amplitude in [("STD", 1e6), ("T1", 1e6), ("T10", 1e4)]:
            for replicate in [1, 2, 3]:
                scans = []
                for i in range(25):
                    rt = 100 + i * 0.5
                    intensity = amplitude * np.exp(-0.5 * ((rt - 106) / 1.0) ** 2)
                    ratio = 0.0107 / 0.9893
                    peaks = [
                        target,
                        intensity,
                        target + 1.003354835,
                        intensity * carbons * ratio,
                        target + 2.00670967,
                        intensity * carbons * (carbons - 1) / 2 * ratio**2,
                    ]
                    encoded = base64.b64encode(struct.pack(">6d", *peaks)).decode("ascii")
                    scans.append(
                        f'<scan num="{i + 1}" msLevel="1" peaksCount="3" polarity="+" '
                        f'retentionTime="PT{rt}S"><peaks precision="64" byteOrder="network" '
                        f'pairOrder="m/z-int">{encoded}</peaks></scan>'
                    )
                xml = (
                    '<?xml version="1.0"?><mzXML xmlns="http://sashimi.sourceforge.net/'
                    'schema_revision/mzXML_3.2"><msRun scanCount="25">'
                    + "".join(scans)
                    + "</msRun></mzXML>"
                )
                path = raw / f"{group}_{pool}_{replicate}.mzXML"
                if path.exists():
                    raise FileExistsError(f"Example file already exists: {path}")
                path.write_text(xml, encoding="utf8")
    library = directory / "library.csv"
    pd.DataFrame(rows).to_csv(library, index=False)
    return raw, library


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    generate_example(args.out)
    print(f"Generated synthetic two-pool example in {args.out}")
