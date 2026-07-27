"""
Developed and tested with Python 3.12.12 and NumPy 2.4.2.

Dresser algorithm for "unspinning" and drift-correcting time-lapse datasets of intranuclear fluorescent spots.

The approach to unspinning is similar to the global Procrustes analyses (GPA) in Python libraries prince and qc-procrustes but is specialized to the nuclear-spot problem by constraining the rotation to be about the measured nuclear center, not immediately solving for translation, not allowing scaling or reflection, and by not using centroids to center the spot clouds (because the spots used for registration generally are off-center; centroids are used to correct for drift).  Correction for nuclear drift is done separately, solved as a slow linear translation of the nuclear center over time, rather than as a per-frame translation.

Rationale
---------
The registration spots - here, specifically 3 Hoechst (H)-emphasized heterochromatin blobs - are adjacent to the inner nuclear surface but are relatively large, can change shape over the course of the time-lapse and their measured coordinate is the center of the blob, not the surface-adjacent position. Thus, each H spot's radial coordinate is noisy and biased inward, while its direction from the nuclear center is reliable. Two ideas drive the logic in this module:

1. The per-frame blob-reshape and measurement noise is averaged out by aligning every frame to a single consensus registration-spot arrangement - the per-spot mean over all aligned frames) - recalculated each iteration. Referencing every frame to the consensus rather than chaining across timepoints reduces the accumulation of calculation drift.

2. Output positions are not projected onto the measurement-determined sphere by default.  Registration, "extra" H and GFP (G) spots marking specific chromosome positions alike keep their real transformed radial position, so the two types of locations are sensibly compared.

Any projection is applied at output, never before the fit.  Rotation fits use the raw registration spot positions about the measured center.  Averaging over all frames in the consensus suppresses radial noise in the fit.

Assumptions from experimental conditions
----------------------------------------
1. The nuclear radius is stable over a timelapse, so one measured "R" serves the whole movie.
2. The nuclear center drifts slowly over the movie.  Spin is removed by a per-frame rotation about the measured center; drift is removed as a single smooth, linear-in-time, translation, not a per-frame one. Estimating a separate per-frame translation from three off-center registration spots is unreliable.  As the nucleus spins, an off-center spot cluster's centroid swings, which a free per-frame translation misreads as drift and injects as spurious motion into the non-registration spot residuals.  A slow linear drift model maintains the real translation and avoids the spurious motion.

Algorithm (for one nucleus)
---------------------------
1. Input the data from a comma-separated values (csv) file with one row per spot per timepoint, giving the spot's label, x,y,z coordinates and timepoint.  Spot labels end in "<spot number>H" or "<spot number>G".  The nuclear center and radius are measured at timepoint 0 (TP0) and given as command-line arguments.  The Hoechst (H) spots used for registration are either specified by their spot numbers or all H spots are used when the cell has exactly three.
2. Subtract the measured "--center" from every spot at every time-point so that the nuclear center sits at the origin which is then the fixed rotation pivot.
3. Pick the three registration H spots, "--H-spots [A,B,C]", or all H spots when the cell has exactly three, and initialize the consensus from their timepoint 0 positions.  Thus, TP0 is the anchor, or reference, orientation.
4. Use a Rotation-only GPA loop until the per-frame rotations stop changing.
   a. For each frame "k" solve the rotation "R_k" about the origin (Kabsch / Wahba singular value decomposition, no centroid removal) best mapping that frame's three registration spots onto the consensus.
   b. Recalculate the consensus as the per-spot mean of every frame's rotated registration spots.
5. Gauge-fix so the TP0 rotation is the identity (output is anchored to TP0).
6. Estimate slow drift by fitting a straight line in time to the spin-removed registration centroid and keep only its time-varying part "drift_k" (0 at TP0).  Each frame's transform is "p -> R_k @ (p - drift_k)", which subtracts for drift then matrix-multiplies by the rotation to calculate the transformed position.  The drift is subtracted from the centroids and then recalculated as the slope of the best-fit line through the spin-removed centroids, so the constant part of the drift (the static offset of the registration cluster from the nuclear center) is kept and only the time-varying part is removed.
7. Apply each transform to all spots (H and G) at frame "k".  The same nuclear rotation and slow drift are removed from the GFP spot, leaving its active motion as the residual.
8. Output the transformed positions without sphere projection so that H and G spots alike keep their real transformed radial position and can be compared without radial position bias. Alternatively, non-H spots can be projected onto the nuclear/sphere surface using the flag "--proj-to-periph", either for spots known to lie on the periphery or to further reduce noise at the expense of exaggerating movements.
9. Output is written to a file with the same name as the input file and partially in the same format as the input but with a "-gpa" suffix, "<stem>-gpa.csv".  Columns are added to the output csv file for each row of spot data, giving the per-frame translation and rotation applied to that spot as "tx,ty,tz" and "rx,ry,rz" respectively.  Rows are added to the output to show the "sphere_radius", the "registration_spots numbers (the H numbers used for registration) and the "track_length' and "avg_speed" per spot.

Command line, flags
----------------------
python -m dresser_unspin_algo "<folder>/<file name>.csv" --center CX CY CZ --radius R

"--center CX CY CZ"  (required) measured nuclear center at timepoint 0.
"--radius R"  (required) measured nuclear radius (assumed stable over the movie).
"--out-dir DIR"  output directory (default: a sister "unspun" folder beside the input, for example "data/<movie name>/unspun/" where the input data are in "data/<movie name>/raw/").
"--time-step-seconds S"  seconds per time-point, used for "avg_speed" (default 60.0; may be given as "60").
"--proj-to-periph"  also project the non-registration (G) spots onto the sphere.  Is off by default so that interior G spots keep their real radial position.
"--H-spots [A,B,C]"  the three H spot numbers to use as registration spots.  Can omit when the cell has exactly three H spots.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import numpy as np

ZERO_TOLERANCE = 1e-9
GPA_MAX_ITERATIONS = 100
GPA_TOLERANCE_RADIANS = 1e-7  # stop when every frame's rotation update is below this

class SpotKind(Enum):
    """Type of intranuclear structure, taken from a label's last letter."""

    HOECHST = "H"
    GFP = "G"

    @classmethod
    def from_label(cls, label: str) -> "SpotKind":
        suffix = label.strip()[-1:].upper()
        for kind in cls:
            if kind.value == suffix:
                return kind
        raise ValueError(
            f"label {label!r} does not end in 'H' or 'G'"
        )


@dataclass(frozen=True)
class Spot:
    """One tracked structure over a contiguous timepoint range.

    "times" is strictly increasing with step 1. "positions" is shaped
    "(len(times), 3)" holding x, y, z in microns, row-aligned to "times".
    """

    label: str
    kind: SpotKind
    times: np.ndarray  # int, shape (T,)
    positions: np.ndarray  # float, shape (T, 3)

    @property
    def start_time(self) -> int:
        return int(self.times[0])

    @property
    def end_time(self) -> int:
        return int(self.times[-1])

    def position_at(self, time: int) -> np.ndarray:
        """Return the (3,) coordinate recorded at "time"."""
        return self.positions[time - self.start_time]


@dataclass(frozen=True)
class Nucleus:
    """All spots measured from a single nucleus, sharing one time range."""

    name: str
    spots: tuple[Spot, ...]
    start_time: int
    end_time: int

    @property
    def time_points(self) -> np.ndarray:
        return np.arange(self.start_time, self.end_time + 1)

    @property
    def num_time_points(self) -> int:
        return self.end_time - self.start_time + 1

    def by_kind(self, kind: SpotKind) -> tuple[Spot, ...]:
        return tuple(s for s in self.spots if s.kind is kind)

    def hoechst_spots(self) -> tuple[Spot, ...]:
        return self.by_kind(SpotKind.HOECHST)

    def gfp_spots(self) -> tuple[Spot, ...]:
        return self.by_kind(SpotKind.GFP)


EXPECTED_HEADER = ["label", "x", "y", "z", "t"]


class ValidationError(Exception):
    """Raised when a CSV departs from the required structure.

    Carries the offending file and (when applicable) spot label so callers can
    report exactly where the data is wrong.
    """

    def __init__(self, path: Path, message: str, label: str | None = None):
        self.path = Path(path)
        self.label = label
        self.detail = message
        where = f"{self.path.name}"
        if label is not None:
            where += f" [spot {label!r}]"
        super().__init__(f"{where}: {message}")


def read_nucleus_csv(path: str | Path) -> Nucleus:
    """Load and validate one nucleus from a CSV file."""
    path = Path(path)
    rows_by_label = _read_rows(path)

    spots: list[Spot] = []
    for label, rows in rows_by_label.items():
        spots.append(_build_spot(path, label, rows))

    _check_shared_time_range(path, spots)

    start = spots[0].start_time
    end = spots[0].end_time
    return Nucleus(
        name=path.stem,
        spots=tuple(spots),
        start_time=start,
        end_time=end,
    )


def _read_rows(path: Path) -> dict[str, list[tuple[int, float, float, float]]]:
    """Parse rows into "{label: [(t, x, y, z), ...]}" in file order."""
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            raise ValidationError(path, "file is empty")

        if [c.strip().lower() for c in header[:5]] != EXPECTED_HEADER:
            raise ValidationError(
                path, f"header must start with {EXPECTED_HEADER}, got {header}"
            )

        rows_by_label: dict[str, list] = defaultdict(list)
        for line_number, row in enumerate(reader, start=2):
            if not row:
                break
            if len(row) < 5:
                raise ValidationError(
                    path, f"line {line_number} has {len(row)} fields, expected >= 5"
                )
            label, raw_x, raw_y, raw_z, raw_t = (c.strip() for c in row[:5])
            try:
                x, y, z = float(raw_x), float(raw_y), float(raw_z)
            except ValueError:
                raise ValidationError(
                    path, f"line {line_number} has non-float coordinates {row[1:4]}",
                    label=label,
                )
            try:
                t = int(raw_t)
            except ValueError:
                raise ValidationError(
                    path, f"line {line_number} has non-integer time-point {raw_t!r}",
                    label=label,
                )
            rows_by_label[label].append((t, x, y, z))

    if not rows_by_label:
        raise ValidationError(path, "file contains a header but no data rows")
    return rows_by_label


def _build_spot(path: Path, label: str, rows: list) -> Spot:
    try:
        kind = SpotKind.from_label(label)
    except ValueError as exc:
        raise ValidationError(path, str(exc), label=label)

    rows_sorted = sorted(rows, key=lambda r: r[0])
    times = [r[0] for r in rows_sorted]

    duplicates = {t for t in times if times.count(t) > 1}
    if duplicates:
        raise ValidationError(
            path, f"duplicate time-point(s) {sorted(duplicates)}", label=label
        )

    expected = list(range(times[0], times[-1] + 1))
    if times != expected:
        missing = sorted(set(expected) - set(times))
        raise ValidationError(
            path,
            f"time-points are not consecutive; missing {missing} "
            f"in range {times[0]}..{times[-1]}",
            label=label,
        )

    positions = np.array([[r[1], r[2], r[3]] for r in rows_sorted], dtype=np.float64)
    return Spot(
        label=label,
        kind=kind,
        times=np.array(times, dtype=np.int64),
        positions=positions,
    )


def _check_shared_time_range(path: Path, spots: list[Spot]) -> None:
    reference = spots[0]
    for spot in spots[1:]:
        if (spot.start_time, spot.end_time) != (
            reference.start_time,
            reference.end_time,
        ):
            raise ValidationError(
                path,
                f"time range {spot.start_time}..{spot.end_time} differs from "
                f"{reference.label!r} range {reference.start_time}..{reference.end_time}",
                label=spot.label,
            )


def _step_rotation(start: np.ndarray, end: np.ndarray) -> np.ndarray | None:
    """Rigid rotation "R" best mapping "start" onto "end" (Kabsch / Wahba).

    "start" and "end" are matched Hoechst point sets "(N, 3)" already
    centered on the sphere center (the rotation center). Returns the proper
    rotation matrix "R" with "R @ start_i ~ end_i", or "None" if the cloud
    has no spatial extent (degenerate; nothing to align).
    """
    if start.shape[0] == 0 or not np.any(np.linalg.norm(start, axis=1) > ZERO_TOLERANCE):
        return None

    covariance = start.T @ end
    u, _, vt = np.linalg.svd(covariance)
    # reflection guard: force a proper rotation (det = +1).
    determinant_sign = np.sign(np.linalg.det(vt.T @ u.T))
    correction = np.diag([1.0, 1.0, determinant_sign])
    return vt.T @ correction @ u.T


def track_stats(
    positions: np.ndarray, time_step_seconds: float
) -> tuple[np.ndarray, np.ndarray]:
    """Per-spot total track length and average speed.

    "positions" is "(num_spots, num_time_points, 3)". "avg_speed" is the
    track length divided by the elapsed time, "(num_time_points - 1) *
    time_step_seconds" (distance units per second).
    """
    steps = np.linalg.norm(np.diff(positions, axis=1), axis=2)  # (spots, steps)
    track_length = steps.sum(axis=1)
    elapsed = max(positions.shape[1] - 1, 1) * time_step_seconds
    return track_length, track_length / elapsed

def _spot_number(label: str) -> int:
    """The spot number embedded in a label, e.g. "cell7Spots 12H" -> "12"."""
    match = re.search(r"(\d+)\s*[A-Za-z]\s*$", label.strip())
    if match is None:
        raise ValueError(f"cannot read a spot number from label {label!r}")
    return int(match.group(1))


def _registration_index(
    nucleus: Nucleus, is_hoechst: list[bool], registration_numbers: list[int] | None
) -> np.ndarray:
    """Spot-stack indices of the H spots used to estimate the rigid transform.

    "None" selects every H spot, which is allowed only when the cell has exactly
    three (more than three needs an explicit "--H-spots" choice). Otherwise the
    requested spot numbers must each name an H spot; a non-H or missing number is
    an error.
    """
    if registration_numbers is None:
        hoechst_index = np.flatnonzero(is_hoechst)
        if hoechst_index.size != 3:
            raise ValueError(
                f"{nucleus.name}: --H-spots omitted but the cell has "
                f"{hoechst_index.size} Hoechst (H) spots, not 3; pass --H-spots "
                f"[A,B,C] to choose the three registration spots"
            )
        return hoechst_index
    number_to_index = {
        _spot_number(s.label): i
        for i, s in enumerate(nucleus.spots)
        if is_hoechst[i]
    }
    missing = [n for n in registration_numbers if n not in number_to_index]
    if missing:
        raise ValueError(
            f"{nucleus.name}: --H-spots {missing} are not Hoechst (H) spot numbers; "
            f"available H spots: {sorted(number_to_index)}"
        )
    return np.array([number_to_index[n] for n in registration_numbers])


def _rotation_vector(rotation: np.ndarray) -> np.ndarray:
    """Axis-angle vector (radians, magnitude = angle) for a 3x3 rotation matrix.
    """
    cos_angle = (np.trace(rotation) - 1.0) / 2.0
    angle = float(np.arccos(np.clip(cos_angle, -1.0, 1.0)))

    if angle < ZERO_TOLERANCE:
        return np.zeros(3)

    axis = np.array(
        [
            rotation[2, 1] - rotation[1, 2],
            rotation[0, 2] - rotation[2, 0],
            rotation[1, 0] - rotation[0, 1],
        ]
    )
    norm = np.linalg.norm(axis)
    if norm < ZERO_TOLERANCE:
        axis = np.sqrt(np.clip((np.diag(rotation) + 1.0) / 2.0, 0.0, None))
        norm = np.linalg.norm(axis)

    return angle * axis / norm


def _project_to_sphere(points: np.ndarray, radius: float) -> np.ndarray:
    """Project "points" radially onto the origin-centered sphere of "radius".
    """

    norms = np.linalg.norm(points, axis=-1, keepdims=True)

    norms = np.where(norms < ZERO_TOLERANCE, 1.0, norms)

    return points / norms * radius


def _gpa_rotations(
    points: np.ndarray, max_iterations: int, tolerance: float
) -> tuple[list[np.ndarray], np.ndarray]:
    """Per-frame **rotations about the measured center** aligning H to a consensus.
    """
    num_time_points = points.shape[1]
    consensus = points[:, 0, :].copy()  # TP0 anchor
    rotations = [np.eye(3)] * num_time_points

    # Generalized Procrustes: alternate two steps until they stop changing.
    #   (a) given the target shape, rotate every frame onto it;
    #   (b) given the rotations, rebuild the target as their average.
    for _ in range(max_iterations):
        new_rotations = [
            _step_rotation(points[:, k, :], consensus) for k in range(num_time_points)
        ]
        new_rotations = [r if r is not None else np.eye(3) for r in new_rotations]
        aligned = np.stack(
            [np.einsum("hc,dc->hd", points[:, k, :], new_rotations[k])
             for k in range(num_time_points)],
            axis=1,
        )
        consensus = aligned.mean(axis=1)

        # Convergence
        max_delta = max(
            np.linalg.norm(_rotation_vector(new @ old.T))
            for new, old in zip(new_rotations, rotations)
        )
        rotations = new_rotations
        if max_delta < tolerance:
            break

    # Gauge-fix: Anchor by making frame 0 the identity
    inverse_0 = rotations[0].T  # rotation matrices are orthogonal: inverse = transpose
    return [inverse_0 @ rotation for rotation in rotations], consensus


def _smooth_drift(
    points: np.ndarray, rotations: list[np.ndarray], consensus: np.ndarray
) -> np.ndarray:
    """Slow (linear-in-time) drift of the nuclear center, in the world frame.
    """
    num_time_points = points.shape[1]
    raw_centroid = points.mean(axis=0)          # (T, 3) registration centroid/frame
    target = consensus.mean(axis=0)             # consensus centroid M
    drift_proxy = np.stack(
        [raw_centroid[k] - rotations[k].T @ target for k in range(num_time_points)]
    )
    times = np.arange(num_time_points)
    design = np.vstack([times, np.ones(num_time_points)]).T
    drift = np.zeros((num_time_points, 3))
    for axis in range(3):
        slope, _intercept = np.linalg.lstsq(design, drift_proxy[:, axis], rcond=None)[0]
        drift[:, axis] = slope * times          # time-varying part only (0 at TP0)
    return drift


def _gpa_transforms(
    points: np.ndarray,
    max_iterations: int = GPA_MAX_ITERATIONS,
    tolerance: float = GPA_TOLERANCE_RADIANS,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Per-frame transforms: rotate about the measured center, remove drift.

    "points" is "(num_hoechst, num_time_points, 3)" of registration H
    positions. Two decoupled steps: solve the rotation about the fixed measured center ("_gpa_rotations"), then estimate a slow linear nuclear-center drift ("_smooth_drift"). Returns one "(R, t)" per time-point where applying "p -> R @ p + t" unspins a spot; "t" encodes the drift removal, "t_k = -R_k @ drift_k" (i.e. "R_k @ (p - drift_k))".
    Subtract world-frame drift first, then undo the spin.
    """
    rotations, consensus = _gpa_rotations(points, max_iterations, tolerance)
    drift = _smooth_drift(points, rotations, consensus)
    return [
        (rotation, -rotation @ drift_k)
        for rotation, drift_k in zip(rotations, drift)
    ]


def gpa_unspin(
    nucleus: Nucleus,
    center: np.ndarray,
    radius: float,
    proj_to_periph: bool = False,
    registration_numbers: list[int] | None = None,
) -> tuple[np.ndarray, list[tuple[int, np.ndarray, np.ndarray]]]:
    """Run the global-Procrustes unspin against a measured center and radius.

    Returns "(positions, transforms)" where "positions" is
    "(num_spots, num_time_points, 3)" of output coordinates and "transforms"
    lists "(time, translation_xyz, rotation_vector_degrees)" per time-point.

    The rigid transform is estimated from a registration subset of three H
    spots. "registration_numbers" selects them by spot number (e.g.
    "[2, 3, 4]"); when "None" every H spot is used (the cell-has-exactly-3-H
    case). Either way the transform is applied to all spots.

    No sphere projection is applied by default: H and G spots alike keep their
    real transformed radial position, so the two locations can be compared directly. "proj_to_periph" projects the G spots onto the measured sphere.
    """
    center = np.asarray(center, dtype=np.float64)

    # Flag which spots are Hoechst (H). Only H spots ride on the nuclear
    # periphery and move rigidly with the nucleus, so only they are trustworthy
    # for measuring the spin. G spots move on their own and would bias the fit.
    is_hoechst = [s.kind is SpotKind.HOECHST for s in nucleus.spots]
    hoechst_index = np.flatnonzero(is_hoechst)
    if hoechst_index.size == 0:
        raise ValueError(f"{nucleus.name}: no Hoechst (H) spots to estimate rotation")
    # Of the H spots, pick the 3 registration spots that actually drive the fit.
    registration_index = _registration_index(nucleus, is_hoechst, registration_numbers)

    # Stack all spots into one array (num_spots, num_time_points, 3), time-aligned.
    positions = np.stack([s.positions for s in nucleus.spots]).astype(np.float64)

    # Translate the whole movie so the measured nuclear center sits at the origin.
    # Every later rotation is about the origin, so this must happen first.
    positions -= center

    # Solve, per frame, the rotation about the measured center plus the slow
    # linear drift, from just the 3 registration H spots (GPA above). Returned as
    # (R_k, t_k) with t_k = -R_k @ drift_k, so p -> R_k @ p + t_k = R_k @ (p - drift_k).
    rigid = _gpa_transforms(positions[registration_index])

    # Apply that same per-frame transform to every spot (H and G). Removing the
    # nucleus's own spin and slow drift from the non-registration and G spots
    # leaves only the real motion as the residual -- the whole point of "unspinning".
    for k, (rotation, translation) in enumerate(rigid):
        positions[:, k, :] = (
            np.einsum("sc,dc->sd", positions[:, k, :], rotation) + translation
        )

    # Repackage the per-frame transforms for logging: absolute time, the drift
    # translation, and the rotation as a degrees axis-angle arrow.
    transforms = [
        (
            nucleus.start_time + k,
            translation.copy(),
            _rotation_vector(rotation) * 180.0 / np.pi,
        )
        for k, (rotation, translation) in enumerate(rigid)
    ]

    # Output coordinates: raw transformed positions without sphere projection by
    # default. H and G alike keep their real transformed radial position, so the
    # two locations are compared directly.
    # --proj-to-periph translates spots onto the nucleus surface/sphere when requested
    output = positions.copy()
    if proj_to_periph:
        gfp_index = np.flatnonzero([s.kind is SpotKind.GFP for s in nucleus.spots])
        output[gfp_index] = _project_to_sphere(positions[gfp_index], radius)
    return output, transforms


def output_path(input_path: str | Path, out_dir: str | Path | None = None) -> Path:
    """Output csv path; defaults to a sister "unspun" folder beside the input.

    For example "data/<movie>/raw/<file>.csv" -> "data/<movie>/unspun/<stem>-gpa.csv".
    """
    if out_dir is None:
        out_dir = Path(input_path).parent.parent / "unspun"
    return Path(out_dir) / f"{Path(input_path).stem}-gpa.csv"


def write_csv(
    path: Path,
    nucleus: Nucleus,
    positions: np.ndarray,
    transforms: list[tuple[int, np.ndarray, np.ndarray]],
    radius: float,
    time_step_seconds: float,
    registration_numbers: list[int],
) -> None:
    """Write coordinates + per-row transforms, then radius, registration spots and
    per-spot track stats. The "registration_spots" footer line records the H
    spot numbers used for the rigid fit so visualization tools can flag them."""
    track_length, avg_speed = track_stats(positions, time_step_seconds)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["label", "x", "y", "z", "t", "tx", "ty", "tz", "rx", "ry", "rz"]
        )
        for i, spot in enumerate(nucleus.spots):
            for k, time in enumerate(range(nucleus.start_time, nucleus.end_time + 1)):
                x, y, z = positions[i, k]
                _, translation, rotvec = transforms[k]
                tx, ty, tz = translation
                rx, ry, rz = rotvec
                writer.writerow(
                    [
                        spot.label,
                        f"{x:.6f}", f"{y:.6f}", f"{z:.6f}", time,
                        f"{tx:.6f}", f"{ty:.6f}", f"{tz:.6f}",
                        f"{rx:.6f}", f"{ry:.6f}", f"{rz:.6f}",
                    ]
                )
        writer.writerow([])
        writer.writerow(["sphere_radius", f"{radius:.6f}"])
        writer.writerow(["registration_spots", *registration_numbers])
        for i, spot in enumerate(nucleus.spots):
            writer.writerow(
                [
                    spot.label,
                    "track_length", f"{track_length[i]:.6f}",
                    "avg_speed", f"{avg_speed[i]:.6f}",
                ]
            )


def _parse_h_spots(text: str) -> list[int]:
    """Parse "--H-spots" like "[2,3,4]" into three distinct ints."""
    parts = [p for p in re.split(r"[,\s]+", text.strip().strip("[]").strip()) if p]
    try:
        numbers = [int(p) for p in parts]
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"--H-spots must be three integers, e.g. [2,3,4]; got {text!r}"
        )
    if len(numbers) != 3:
        raise argparse.ArgumentTypeError(
            f"--H-spots needs exactly 3 spot numbers, e.g. [2,3,4]; got {numbers}"
        )
    if len(set(numbers)) != 3:
        raise argparse.ArgumentTypeError(f"--H-spots numbers must be distinct; got {numbers}")
    return numbers


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Unspin via global Procrustes against a measured center and radius."
    )
    parser.add_argument("input", help="raw spot CSV file")
    parser.add_argument(
        "--center", type=float, nargs=3, metavar=("CX", "CY", "CZ"), required=True,
        help="measured nuclear center at timepoint 0 (microns)",
    )
    parser.add_argument(
        "--radius", type=float, required=True,
        help="measured nuclear radius (assumed stable over the movie)",
    )
    parser.add_argument(
        "--out-dir", default=None,
        help="output directory (default: sister 'unspun' folder beside the input)",
    )
    parser.add_argument(
        "--time-step-seconds", type=float, default=60.0,
        help="seconds per timepoint, used for avg_speed (default 60.0)",
    )
    parser.add_argument(
        "--proj-to-periph", dest="proj_to_periph", action="store_true",
        help="project G spots onto the sphere; off by default so "
        "interior G spots keep their real radial position",
    )
    parser.add_argument(
        "--H-spots", dest="h_spots", type=_parse_h_spots, default=None,
        metavar="[A,B,C]",
        help="the 3 H spot numbers to use as registration spots, e.g. [2,3,4]; "
        "can omit when the cell has exactly 3 H spots",
    )
    args = parser.parse_args(argv)

    try:
        nucleus = read_nucleus_csv(args.input)
        positions, transforms = gpa_unspin(
            nucleus, np.array(args.center), args.radius,
            proj_to_periph=args.proj_to_periph,
            registration_numbers=args.h_spots,
        )
    except (ValidationError, ValueError) as exc:
        print(f"ERROR {Path(args.input).name}: {exc}", file=sys.stderr)
        return 1

    # The registration H numbers actually used (all 3 H spots when --H-spots omitted).
    registration_numbers = args.h_spots
    if registration_numbers is None:
        registration_numbers = [
            _spot_number(s.label)
            for s in nucleus.spots
            if s.kind is SpotKind.HOECHST
        ]

    destination = output_path(args.input, args.out_dir)
    write_csv(
        destination, nucleus, positions, transforms, args.radius,
        args.time_step_seconds, registration_numbers,
    )

    cx, cy, cz = args.center
    print(
        f"{Path(args.input).name}: center=({cx:.3f}, {cy:.3f}, {cz:.3f}) "
        f"radius={args.radius:.3f}"
    )
    track_length, avg_speed = track_stats(positions, args.time_step_seconds)
    for spot, length, speed in zip(nucleus.spots, track_length, avg_speed):
        print(f"  {spot.label:>3}  track_length={length:8.3f}  avg_speed={speed:7.3f}")
    print(f"-> {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
