#!/usr/bin/env python

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import types
from pathlib import Path

import nibabel as nb
import numpy as np


def _write_error(out_dir: Path, message: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "error.txt").write_text(message)


def _probe_nifreeze() -> tuple[bool, str]:
    if shutil.which("antsRegistration") is None:
      return False, "Missing ANTs executable: antsRegistration"

    try:
        import pkgutil

        missing = [name for name in ("nitransforms", "nipype") if pkgutil.find_loader(name) is None]
        if missing:
            return False, "Missing Python dependency/dependencies: " + ", ".join(missing)

        vendor_src = Path("vendor/nifreeze/src").resolve()
        sys.path.insert(0, str(vendor_src))
        version_mod = types.ModuleType("nifreeze._version")
        version_mod.__version__ = "0"
        sys.modules["nifreeze._version"] = version_mod

        from nifreeze.data import load  # noqa: F401
        from nifreeze.estimator import Estimator  # noqa: F401
    except Exception as exc:
        return False, f"Failed to import nifreeze stack: {exc!r}"

    return True, "ok"


def _write_motion_csv(mats: np.ndarray, path: Path) -> None:
    fieldnames = ["frame"] + [f"m{i}{j}" for i in range(1, 5) for j in range(1, 5)]
    with path.open("w", newline="") as fobj:
        writer = csv.DictWriter(fobj, fieldnames=fieldnames)
        writer.writeheader()
        for i in range(mats.shape[0]):
            row = {"frame": i + 1}
            values = mats[i].reshape(-1, order="F")
            row.update({fieldnames[j + 1]: float(values[j]) for j in range(16)})
            writer.writerow(row)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--input")
    parser.add_argument("--mask")
    parser.add_argument("--out-dir")
    parser.add_argument("--model", default="mean")
    parser.add_argument("--nthreads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=4401)
    args = parser.parse_args(argv)

    ready, message = _probe_nifreeze()
    if args.probe:
        print(message)
        return 0 if ready else 2

    if not ready:
        _write_error(Path(args.out_dir), message)
        return 2

    vendor_src = Path("vendor/nifreeze/src").resolve()
    sys.path.insert(0, str(vendor_src))
    version_mod = types.ModuleType("nifreeze._version")
    version_mod.__version__ = "0"
    sys.modules["nifreeze._version"] = version_mod

    from nifreeze.data import load
    from nifreeze.estimator import Estimator

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dataset = load(args.input, brainmask_file=args.mask)
    estimator = Estimator(args.model, single_fit=False)
    estimator.run(dataset, omp_nthreads=args.nthreads, n_jobs=1, seed=args.seed)

    corrected_path = out_dir / "corrected.nii.gz"
    dataset.to_nifti(corrected_path, write_hmxfms=False)

    motion_affines = dataset.motion_affines
    if motion_affines is None:
        motion_affines = np.repeat(np.eye(4)[None, ...], len(dataset), axis=0)

    _write_motion_csv(motion_affines, out_dir / "motion_mats.csv")
    with (out_dir / "run_metadata.csv").open("w", newline="") as fobj:
        writer = csv.DictWriter(fobj, fieldnames=["method", "model"])
        writer.writeheader()
        writer.writerow({"method": "nifreeze", "model": args.model})

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
