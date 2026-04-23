#!/usr/bin/env python

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import tempfile
import types
from collections import namedtuple
from pathlib import Path

import nibabel as nb
import numpy as np


def _write_error(out_dir: Path, message: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "error.txt").write_text(message)


def _bootstrap_nifreeze_registration():
    vendor_src = Path("vendor/nifreeze/src").resolve()
    if str(vendor_src) not in sys.path:
        sys.path.insert(0, str(vendor_src))

    version_mod = types.ModuleType("nifreeze._version")
    version_mod.__version__ = "0"
    sys.modules["nifreeze._version"] = version_mod

    from nifreeze.registration.ants import _prepare_registration_data, _run_registration
    from nitransforms.linear import LinearTransformsMapping
    from nitransforms.resampling import apply

    return _prepare_registration_data, _run_registration, LinearTransformsMapping, apply


def _probe_nifreeze() -> tuple[bool, str]:
    if shutil.which("antsRegistration") is None:
        return False, "Missing ANTs executable: antsRegistration"

    try:
        import pkgutil

        missing = [name for name in ("nitransforms", "nipype") if pkgutil.find_loader(name) is None]
        if missing:
            return False, "Missing Python dependency/dependencies: " + ", ".join(missing)

        _bootstrap_nifreeze_registration()
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


def _model_stat(model: str) -> str:
    key = model.lower()
    if key in {"avg", "average", "mean"}:
        # This mirrors nifreeze's current ExpectationModel default.
        return "median"
    raise ValueError(f"Unsupported nifreeze model for this benchmark runner: {model}")


def _expectation_prediction(data: np.ndarray, index: int, stat: str) -> np.ndarray:
    index_mask = np.ones(data.shape[-1], dtype=bool)
    index_mask[index] = False
    return getattr(np, stat)(data[..., index_mask], axis=-1)


def _resample_corrected(
    data: np.ndarray,
    affine: np.ndarray,
    header: nb.Nifti1Header,
    motion_affines: np.ndarray,
):
    _, _, LinearTransformsMapping, apply = _bootstrap_nifreeze_registration()

    image_grid = namedtuple("ImageGrid", ("shape", "affine"))
    reference = image_grid(shape=data.shape[:3], affine=affine)
    xforms = LinearTransformsMapping(motion_affines, reference=reference)
    corrected = np.empty_like(data)

    for i, xform in enumerate(xforms):
        moving = nb.Nifti1Image(data[..., i], affine, header)
        corrected[..., i] = np.asanyarray(apply(xform, moving, order=3).dataobj, dtype=data.dtype)

    return corrected


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

    try:
        stat = _model_stat(args.model)
    except ValueError as exc:
        _write_error(Path(args.out_dir), str(exc))
        return 2

    _prepare_registration_data, _run_registration, _, _ = _bootstrap_nifreeze_registration()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    img = nb.load(args.input)
    data = np.asanyarray(img.dataobj)
    affine = img.affine
    header = img.header.copy()
    motion_affines = np.repeat(np.eye(4)[None, ...], data.shape[-1], axis=0)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        for i in range(data.shape[-1]):
            predicted = _expectation_prediction(data, i, stat=stat)
            predicted_path, volume_path, _ = _prepare_registration_data(
                data[..., i],
                predicted,
                affine,
                i,
                tmp_path,
                "both",
            )
            xform = _run_registration(
                predicted_path,
                volume_path,
                i,
                tmp_path,
                fixedmask_path=args.mask,
                output_transform_prefix=f"ants-{i:05d}",
                num_threads=args.nthreads,
                seed=args.seed,
            )
            motion_affines[i] = xform.matrix

    corrected = _resample_corrected(data, affine, header, motion_affines)
    corrected_path = out_dir / "corrected.nii.gz"
    nb.Nifti1Image(corrected, affine, header).to_filename(corrected_path)

    _write_motion_csv(motion_affines, out_dir / "motion_mats.csv")
    with (out_dir / "run_metadata.csv").open("w", newline="") as fobj:
        writer = csv.DictWriter(fobj, fieldnames=["method", "model", "stat"])
        writer.writeheader()
        writer.writerow({"method": "nifreeze", "model": args.model, "stat": stat})

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
