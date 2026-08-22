#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_log.py — turn galactic_plane_log.csv into science plots.

Companion to hydrogen_mapper.py.  Reads the CSV the mapper writes (plus the
per-scan spectra saved in spectra/*.npz) and produces, next to the CSV:

  1. lv_diagram.png      — longitude–velocity peak scatter with the tangent
                           envelope (each lane of points ≈ a spiral arm).
  2. lv_diagram_2d.png   — TRUE l–v intensity image rebuilt from the saved
                           spectra (every velocity channel, not just peaks).
                           Skipped if no spectra have been saved yet.
  3. rotation_curve.png  — tangent-point rotation curve V(R) from the
                           inner-galaxy sightlines (0<l<90 or 270<l<360).
  4. faceon_map.png      — top-down map of the Milky Way: each (l, v_LSR)
                           converted to a kinematic distance and plotted in
                           galactocentric coordinates (Sun at (0, 8.2) kpc).
                           Inner-galaxy points have a near/far ambiguity, so
                           both solutions are shown (near = filled, far = open).

Rows recorded by the always-record fallback with no real detection are junk;
they are excluded with an SNR cut (default 5 sigma; rows from logs older than
the SNR_sigma column always pass).

Usage
-----
    python analyze_log.py                    # galactic_plane_log.csv next to this file
    python analyze_log.py --show             # also pop the plots up on screen
    python analyze_log.py --hardware-only    # ignore rows tagged Source=Simulation
    python analyze_log.py --min-snr 5        # detection threshold (sigma)
    python analyze_log.py --fits             # also export the log as a FITS table
    python analyze_log.py path/to/log.csv    # a specific log file
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys

import numpy as np

import matplotlib
if "--show" not in sys.argv:            # no display needed when only saving
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Constants (keep consistent with hydrogen_mapper.py).
F_REST_HZ = 1420.405752e6
C_KMS = 299792.458
R0_KPC = 8.2
V0_KMS = 220.0


def load_rows(path: str, hardware_only: bool):
    """Arrays (l, v_lsr, power, snr, fwhm) plus the set of saved spectrum paths.

    snr/fwhm are NaN for rows written before those columns existed.
    """
    l, v, p, s, w = [], [], [], [], []
    spec_paths = set()
    base_dir = os.path.dirname(os.path.abspath(path))
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if hardware_only and row.get("Source", "").strip() == "Simulation":
                continue
            try:
                l.append(float(row["Galactic_Longitude_L"]))
                v.append(float(row["Doppler_Velocity_LSR_km_s"]))
                p.append(float(row["Peak_Power_dB"]))
            except (KeyError, ValueError):
                continue
            try:
                s.append(float(row.get("SNR_sigma", "")))
            except (TypeError, ValueError):
                s.append(np.nan)
            try:
                w.append(float(row.get("Line_FWHM_kHz", "")))
            except (TypeError, ValueError):
                w.append(np.nan)
            sf = (row.get("Spectrum_File") or "").strip()
            if sf and sf != "unsaved":
                full = os.path.join(base_dir, sf)
                if os.path.exists(full):
                    spec_paths.add(full)
    return (np.array(l), np.array(v), np.array(p), np.array(s), np.array(w),
            sorted(spec_paths))


# Any real HI feature is at least a few km/s wide (>=~20 kHz).  A recorded
# width below this means the always-record fallback landed on a noise spike or
# an unresolved bandpass artifact, not a line.
MIN_FWHM_KHZ = 20.0


def snr_mask(snr: np.ndarray, min_snr: float,
             fwhm: np.ndarray | None = None) -> np.ndarray:
    """Keep rows at/above the cuts; legacy rows without the columns pass."""
    keep = np.isnan(snr) | (snr >= min_snr)
    if fwhm is not None:
        keep &= np.isnan(fwhm) | (fwhm >= MIN_FWHM_KHZ)
    return keep


def tangent_envelope(ax):
    """Overlay the flat-rotation tangent-velocity envelope on an l–v axes."""
    lg = np.linspace(0.5, 89.5, 200)
    ax.plot(lg, V0_KMS * (1.0 - np.sin(np.radians(lg))), ls="--", lw=0.9,
            color="#888888", label="tangent envelope (flat curve)")
    lg4 = np.linspace(270.5, 359.5, 200)
    ax.plot(lg4, -V0_KMS * (1.0 - np.abs(np.sin(np.radians(lg4)))),
            ls="--", lw=0.9, color="#888888")


def plot_lv(l, v, p, out_path: str):
    """Longitude–velocity diagram: ridges/lanes of points trace spiral arms."""
    fig, ax = plt.subplots(figsize=(9, 5))
    sc = ax.scatter(l, v, c=p, cmap="plasma", s=40, edgecolors="k", linewidths=0.3)
    ax.axhline(0, color="gray", lw=0.6)
    tangent_envelope(ax)
    ax.set_xlim(0, 360)
    ax.set_ylim(-250, 250)
    ax.set_xlabel("Galactic Longitude  l  (deg)")
    ax.set_ylabel("v_LSR  (km/s)")
    ax.set_title("Longitude–Velocity Diagram  (each lane ≈ a spiral arm)")
    ax.legend(loc="upper right", fontsize=8)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Peak Power (dB)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote {out_path}")
    return fig


def plot_lv_image(spec_paths, l_pts, v_pts, out_path: str):
    """TRUE l–v diagram: a 2-D intensity image built from the saved spectra.

    Each scan contributes a full vertical stripe of velocity channels at its
    longitude — arms appear as continuous bright lanes, including blended
    components the peak-finder would miss.  Returns the figure, or None if no
    spectra are available.
    """
    v_lo, v_hi, dv = -250.0, 250.0, 2.0
    nv = int((v_hi - v_lo) / dv)
    img = np.full((nv, 360), np.nan)
    used = 0
    for path in spec_paths:
        try:
            with np.load(path) as d:
                f = d["freqs_hz"]
                psd = d["psd_lin"]
                base = d["baseline_lin"] if "baseline_lin" in d.files else np.array([])
                l_deg = float(d["l"])
                lsr = float(d["lsr_corr"])
        except Exception as exc:
            print(f"  skipping unreadable spectrum {os.path.basename(path)}: {exc}")
            continue
        if base.size == psd.size and base.size > 0:
            t = (psd - base) / np.maximum(base, 1e-30)   # line-to-continuum
        else:
            t = psd / np.median(psd) - 1.0               # crude flattening
        v = C_KMS * (F_REST_HZ - f) / F_REST_HZ + lsr
        idx = ((v - v_lo) / dv).astype(int)
        ok = (idx >= 0) & (idx < nv)
        col = np.full(nv, np.nan)
        np.fmax.at(col, idx[ok], t[ok])
        li = int(l_deg) % 360
        img[:, li] = np.fmax(img[:, li], col)
        used += 1
    if used == 0:
        return None

    fig, ax = plt.subplots(figsize=(9, 5))
    finite = img[np.isfinite(img)]
    vmax = float(np.percentile(finite, 99.5)) if finite.size else 1.0
    masked = np.ma.masked_invalid(img)
    cmap = plt.get_cmap("inferno").copy()
    cmap.set_bad("#101010")
    im = ax.imshow(masked, origin="lower", aspect="auto",
                   extent=[0, 360, v_lo, v_hi], cmap=cmap,
                   vmin=0.0, vmax=max(vmax, 1e-3), interpolation="nearest")
    if len(l_pts):
        ax.scatter(l_pts, v_pts, s=14, marker="x", color="#66ccff",
                   linewidths=0.8, label="logged peaks")
        ax.legend(loc="upper right", fontsize=8)
    ax.set_xlabel("Galactic Longitude  l  (deg)")
    ax.set_ylabel("v_LSR  (km/s)")
    ax.set_title(f"l–v Intensity Diagram  ({used} scan(s), full line profiles)")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("(ON − OFF) / OFF  line-to-continuum ratio")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote {out_path}")
    return fig


def plot_rotation_curve(l, v, out_path: str):
    """Tangent-point rotation curve from inner-galaxy sightlines.

    For a longitude l in the inner galaxy the most extreme LSR velocity comes
    from the tangent point, where the galactocentric radius is R = R0*|sin l|
    and the rotation speed is  V(R) = |v_tangent| + V0*|sin l|.
    """
    # Bin by integer longitude and take the extreme (tangent) velocity per bin.
    Rs, Vs = [], []
    for lon in np.unique(np.round(l).astype(int)):
        sel = np.round(l).astype(int) == lon
        vv = v[sel]
        if 0 < lon < 90:               # Quadrant I: tangent = max (reddest)
            v_tan = np.max(vv)
        elif 270 < lon < 360:          # Quadrant IV: tangent = min (bluest)
            v_tan = np.min(vv)
        else:
            continue                   # outer galaxy: tangent method invalid
        sinl = abs(np.sin(np.radians(lon)))
        if sinl < 0.05:
            continue                   # near l=0/180: degenerate
        Rs.append(R0_KPC * sinl)
        Vs.append(abs(v_tan) + V0_KMS * sinl)

    has_points = len(Rs) > 0
    fig, ax = plt.subplots(figsize=(8, 5))
    if has_points:
        order = np.argsort(Rs)
        Rs, Vs = np.array(Rs)[order], np.array(Vs)[order]
        ax.plot(Rs, Vs, "o-", color="#1f77b4", label="measured (tangent point)")
    ax.axhline(V0_KMS, color="gray", ls="--", lw=0.8,
               label=f"V0 = {V0_KMS:.0f} km/s")
    ax.set_xlabel("Galactocentric radius  R  (kpc)")
    ax.set_ylabel("Rotation speed  V(R)  (km/s)")
    ax.set_title("Galactic Rotation Curve (tangent-point method)")
    ax.set_xlim(0, R0_KPC)
    ax.set_ylim(0, 300)                 # fixed, readable scale (no offset notation)
    ax.yaxis.get_major_formatter().set_useOffset(False)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote {out_path}")
    if not has_points:
        print("  (no inner-galaxy points yet — collect scans at 0°<l<90° "
              "or 270°<l<360° to build the rotation curve)")
    return fig


def kinematic_solutions(l_deg: float, v_kms: float):
    """(x, y, kind) galactocentric positions for one (l, v_LSR) measurement.

    Flat rotation curve:  v = V0 * sin(l) * (R0/R - 1)
      ->  R = R0 / (1 + v / (V0 sin l))
    Distance along the sightline:  r = R0 cos l ± sqrt(R² − R0² sin²l).
    Inside the solar circle both roots are physical (near/far ambiguity);
    outside, only one root is positive.  Sun at (0, R0), GC at the origin.
    """
    lr = math.radians(l_deg)
    sinl, cosl = math.sin(lr), math.cos(lr)
    if abs(sinl) < 0.05:                       # l ~ 0/180: v ~ 0, degenerate
        return []
    denom = 1.0 + v_kms / (V0_KMS * sinl)
    if denom <= 0.05:                          # unphysical (forbidden velocity)
        return []
    R = R0_KPC / denom
    if not (0.5 <= R <= 25.0):
        return []
    disc = R * R - (R0_KPC * sinl) ** 2
    root = math.sqrt(disc) if disc > 0 else 0.0
    out = []
    tangent = root < 0.05
    for r, kind in ((R0_KPC * cosl - root, "near"), (R0_KPC * cosl + root, "far")):
        if r <= 0.05:
            continue
        if tangent:
            kind = "tangent"
        x = r * sinl
        y = R0_KPC - r * cosl
        out.append((x, y, kind))
        if tangent:                            # both roots identical: keep one
            break
    if len(out) == 1 and out[0][2] == "far" and R > R0_KPC:
        out[0] = (out[0][0], out[0][1], "outer")   # unique outer-galaxy solution
    return out


def plot_faceon(l, v, out_path: str):
    """Top-down spiral-arm map from kinematic distances."""
    pts = {"near": [], "far": [], "tangent": [], "outer": []}
    vs = {"near": [], "far": [], "tangent": [], "outer": []}
    for li, vi in zip(l, v):
        for x, y, kind in kinematic_solutions(float(li), float(vi)):
            pts[kind].append((x, y))
            vs[kind].append(vi)

    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    # Galactic scaffolding: GC, Sun, solar circle, radius rings.
    for rad in (4, 12):
        ax.add_patch(plt.Circle((0, 0), rad, fill=False, color="#cccccc",
                                lw=0.6, ls=":"))
    ax.add_patch(plt.Circle((0, 0), R0_KPC, fill=False, color="#999999",
                            lw=0.8, ls="--"))
    ax.plot(0, 0, "*", ms=16, color="#e6b400", label="Galactic Center")
    ax.plot(0, R0_KPC, "o", ms=9, mfc="#4da3ff", mec="k", label="Sun")

    style = {"near":    dict(marker="o", filled=True,  label="near solution"),
             "far":     dict(marker="o", filled=False, label="far solution"),
             "tangent": dict(marker="s", filled=True,  label="tangent point"),
             "outer":   dict(marker="D", filled=True,  label="outer galaxy")}
    cmap = plt.get_cmap("coolwarm")
    norm = matplotlib.colors.Normalize(vmin=-150, vmax=150)
    sc = None
    for kind, st in style.items():
        if not pts[kind]:
            continue
        xy = np.array(pts[kind])
        kw = dict(marker=st["marker"], s=45, zorder=3, label=st["label"])
        if st["filled"]:
            sc = ax.scatter(xy[:, 0], xy[:, 1], c=vs[kind], cmap=cmap,
                            norm=norm, edgecolors="k", linewidths=0.4, **kw)
        else:
            # open markers: velocity colour on the EDGE (c= would fill them)
            ax.scatter(xy[:, 0], xy[:, 1], facecolors="none",
                       edgecolors=cmap(norm(np.asarray(vs[kind]))),
                       linewidths=1.2, **kw)
    if sc is None:                       # only open markers plotted: still need
        sc = ax.scatter([], [], c=[], cmap=cmap, norm=norm)  # a colorbar handle
    ax.set_xlabel("X  (kpc)")
    ax.set_ylabel("Y  (kpc)   —   GC at origin, Sun at (0, 8.2)")
    ax.set_title("Face-on Milky Way Map (kinematic distances, flat rotation)")
    ax.set_aspect("equal")
    ax.set_xlim(-14, 14)
    ax.set_ylim(-6, 16)
    ax.legend(loc="lower left", fontsize=8)
    if sc is not None:
        cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
        cb.set_label("v_LSR (km/s)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote {out_path}")
    n_pts = sum(len(x) for x in pts.values())
    if n_pts == 0:
        print("  (no plottable kinematic points yet)")
    return fig


def export_fits(csv_path: str, out_path: str):
    """Optional FITS-table export of the whole log (astronomy interchange)."""
    from astropy.table import Table
    t = Table.read(csv_path, format="ascii.csv")
    t.write(out_path, format="fits", overwrite=True)
    print(f"wrote {out_path}")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    default_csv = os.path.join(here, "galactic_plane_log.csv")

    ap = argparse.ArgumentParser(description="Plot the HI drift-scan log.")
    ap.add_argument("csv", nargs="?", default=default_csv,
                    help="path to galactic_plane_log.csv")
    ap.add_argument("--show", action="store_true", help="display plots on screen")
    ap.add_argument("--hardware-only", action="store_true",
                    help="ignore rows tagged Source=Simulation")
    ap.add_argument("--min-snr", type=float, default=5.0,
                    help="drop rows below this detection SNR in sigma "
                         "(legacy rows without an SNR always pass; default 5)")
    ap.add_argument("--fits", action="store_true",
                    help="also export the log as a FITS table")
    args = ap.parse_args()

    if not os.path.exists(args.csv):
        print(f"No log file at {args.csv}")
        sys.exit(1)

    l, v, p, snr, fwhm, spec_paths = load_rows(args.csv, args.hardware_only)
    if len(l) == 0:
        print("No usable rows in the log yet.")
        sys.exit(1)

    keep = snr_mask(snr, args.min_snr, fwhm)
    dropped = int(np.sum(~keep))
    l, v, p = l[keep], v[keep], p[keep]
    print(f"Loaded {len(l)} measurement(s) from {os.path.basename(args.csv)}"
          + (f" ({dropped} below {args.min_snr:.0f}-sigma or narrower than "
             f"{MIN_FWHM_KHZ:.0f} kHz dropped)" if dropped else ""))
    if len(l) == 0:
        print("All rows fell below the SNR cut — nothing to plot.")
        sys.exit(1)

    out_dir = os.path.dirname(os.path.abspath(args.csv))
    plot_lv(l, v, p, os.path.join(out_dir, "lv_diagram.png"))
    if plot_lv_image(spec_paths, l, v,
                     os.path.join(out_dir, "lv_diagram_2d.png")) is None:
        print("no saved spectra yet — skipping lv_diagram_2d.png "
              "(new scans save them automatically)")
    plot_rotation_curve(l, v, os.path.join(out_dir, "rotation_curve.png"))
    plot_faceon(l, v, os.path.join(out_dir, "faceon_map.png"))

    if args.fits:
        export_fits(args.csv, os.path.join(out_dir, "galactic_plane_log.fits"))

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
