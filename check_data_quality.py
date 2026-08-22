#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data-quality checks on the saved drift-scan spectra.

Everything reported in the README's "Is the signal real?" section is produced
here, from the committed `spectra/*.npz` alone.  No hardware, no network:

    python check_data_quality.py

The four tests, in order of what they are good for:

1. SKY-LOCK -- the one that matters.  A real astronomical line sits at a fixed
   velocity in the Local Standard of Rest and therefore at a *moving* radio
   frequency, because the observer's own velocity changes through the night and
   through the year.  An instrumental artefact does the opposite: fixed in
   frequency, drifting in LSR velocity.  Measuring the peak in both frames and
   correlating each against the applied LSR correction separates the two
   without needing any absolute calibration.

2. VELOCITY EXTENT -- how far out in velocity emission is actually detected,
   against the tangent-point velocity expected for a flat rotation curve.  This
   is what decides whether a tangent-point rotation curve is a measurement or
   just `V0*sin(l)` rearranged.

3. BASELINE DRIFT -- the off-line continuum ratio, which is zero when the
   cold-sky reference is still valid and grows as the receiver warms.

4. DRIFT vs FLUX -- guards against the obvious self-deception in test 3: if the
   apparent line were being manufactured by that drift, line flux would climb
   with it.
"""

from __future__ import annotations

import glob
import os
import sys

import numpy as np

F_REST_HZ = 1420.405752e6
C_KMS = 299792.458
V0_KMS = 220.0                  # flat-rotation circular speed
SMOOTH_BINS = 85                # ~50 kHz, matched to the HI line width
LINE_WINDOW = (-40.0, 60.0)     # km/s, where the bright local component sits
QUIET_WINDOW = (200.0, 300.0)   # km/s, no galactic HI expected out here


def load(path: str):
    """Return (v_lsr, v_topo, ratio, l, b, integration_s) sorted by velocity."""
    d = np.load(path, allow_pickle=True)
    freqs = np.asarray(d["freqs_hz"])
    psd = np.asarray(d["psd_lin"])
    base = np.asarray(d["baseline_lin"])

    # Line-to-continuum ratio.  Dividing by the OFF scan removes the receiver
    # bandpass; what survives is (sky_on - sky_off) / sky_off.
    ratio = (psd - base) / base
    v_topo = C_KMS * (F_REST_HZ - freqs) / F_REST_HZ
    v_lsr = v_topo + float(d["lsr_corr"])

    order = np.argsort(v_lsr)
    return (v_lsr[order], v_topo[order], ratio[order],
            float(d["l"]), float(d["b"]), float(d["integration_s"]),
            float(d["lsr_corr"]))


def smooth(y: np.ndarray, n: int = SMOOTH_BINS) -> np.ndarray:
    return np.convolve(y, np.ones(n) / n, mode="same")


def detrend(v: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Remove the linear continuum tilt left by gain drift.

    Without this, a scan whose continuum sits 40 % low swamps any line
    measurement; the tilt is instrumental and carries no velocity information.
    """
    w = (v > -130) & (v < 300)
    return y - np.polyval(np.polyfit(v[w], y[w], 1), v)


def robust_rms(y: np.ndarray) -> float:
    """MAD-based noise estimate: immune to the line itself and to spikes."""
    return float(1.4826 * np.median(np.abs(y - np.median(y))))


def _trapz(y: np.ndarray, x: np.ndarray) -> float:
    """np.trapz was renamed np.trapezoid in numpy 2.0 and removed in 2.3."""
    fn = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return float(fn(y, x))


def terminal_velocity(v: np.ndarray, rs: np.ndarray, rms: float) -> float:
    """Upper velocity edge of the emission, measured properly.

    Take the contiguous run of channels above 3 sigma that contains the peak
    and return its high-velocity end.  Using the highest 3-sigma channel
    anywhere in the window instead would latch onto isolated noise spikes far
    from the line and badly overstate the extent.
    """
    m = (v > -20) & (v < 260)
    idx = np.flatnonzero(m)
    if idx.size == 0:
        return float("nan")
    peak = idx[int(np.argmax(rs[idx]))]
    if rs[peak] < 3 * rms:
        return float("nan")
    i = peak
    while i + 1 < rs.size and rs[i + 1] > 3 * rms:
        i += 1
    return float(v[i])


def test_sky_lock(scans) -> None:
    print("\n" + "=" * 78)
    print(" TEST 1  Is the line locked to the sky or to the receiver?")
    print("=" * 78)
    print(" A real line holds still in the LSR frame while its radio frequency moves.")
    print(" An artefact holds still in frequency while its LSR velocity moves.\n")
    print(f" {'l':>6} {'b':>6} {'lsr_corr':>9} {'peak v_LSR':>11} {'peak v_topo':>12}")

    lsr, p_lsr, p_topo = [], [], []
    for v, vt, r, l, b, _it, corr in scans:
        rs = detrend(v, smooth(r))
        m = (v > -60) & (v < 60)
        i = int(np.argmax(np.where(m, rs, -np.inf)))
        print(f" {l:6.1f} {b:+6.1f} {corr:+9.1f} {v[i]:+11.1f} {vt[i]:+12.1f}")
        lsr.append(corr)
        p_lsr.append(v[i])
        p_topo.append(vt[i])

    lsr, p_lsr, p_topo = map(np.asarray, (lsr, p_lsr, p_topo))
    c_lsr = float(np.corrcoef(lsr, p_lsr)[0, 1])
    c_topo = float(np.corrcoef(lsr, p_topo)[0, 1])

    print(f"\n LSR correction spanned {lsr.max() - lsr.min():.1f} km/s across the data set.")
    print(f"   peak in LSR (sky) frame      : {p_lsr.mean():+6.1f} +/- {p_lsr.std():.1f} km/s")
    print(f"   peak in topocentric frame    : {p_topo.mean():+6.1f} +/- {p_topo.std():.1f} km/s")
    print(f"   corr(lsr_corr, peak v_LSR)   : {c_lsr:+.2f}   (0 = sky-locked)")
    print(f"   corr(lsr_corr, peak v_topo)  : {c_topo:+.2f}   (0 = frequency-locked)")
    verdict = ("SKY-LOCKED -- consistent with real galactic HI"
               if abs(c_lsr) < abs(c_topo) else
               "FREQUENCY-LOCKED -- consistent with an instrumental artefact")
    print(f"\n   VERDICT: {verdict}")


def test_velocity_extent(scans) -> None:
    print("\n" + "=" * 78)
    print(" TEST 2  How far out in velocity is emission actually detected?")
    print("=" * 78)
    print(" v_tangent = V0*(1 - sin l) is where a flat rotation curve puts the")
    print(" terminal velocity.  Falling short means the tangent point was never")
    print(" measured, and a tangent-point rotation curve would be circular.\n")
    print(f" {'l':>6} {'integ':>6} {'v_tan expected':>15} {'v extent measured':>18}  ")

    for v, _vt, r, l, _b, it, _c in scans:
        rs = detrend(v, smooth(r))
        rms = robust_rms(rs[(v > QUIET_WINDOW[0]) & (v < QUIET_WINDOW[1])])
        v_tan = V0_KMS * (1.0 - np.sin(np.radians(l)))
        extent = terminal_velocity(v, rs, rms)
        mark = ("no detection" if not np.isfinite(extent)
                else "reaches it" if extent >= v_tan else "short")
        print(f" {l:6.1f} {it:6.0f} {v_tan:+15.1f} {extent:+18.1f}  {mark}")


def test_baseline_drift(scans) -> list[float]:
    print("\n" + "=" * 78)
    print(" TEST 3  How stale was the cold-sky baseline?")
    print("=" * 78)
    print(" Off-line continuum ratio.  0 = the OFF reference still describes the")
    print(" receiver.  Large negative = gain drifted away since it was captured.\n")
    print(f" {'l':>6} {'integ':>6} {'continuum':>10}  verdict")

    drifts = []
    for v, _vt, r, l, _b, it, _c in scans:
        cont = float(np.median(r[(v > 120) & (v < 300)]))
        verdict = ("ok" if abs(cont) < 0.05
                   else "drifting" if abs(cont) < 0.25 else "SEVERE")
        print(f" {l:6.1f} {it:6.0f} {cont:+10.3f}  {verdict}")
        drifts.append(cont)
    return drifts


def test_drift_vs_flux(scans, drifts: list[float]) -> None:
    print("\n" + "=" * 78)
    print(" TEST 4  Is the drift in test 3 manufacturing the line?")
    print("=" * 78)
    print(" If it were, integrated line flux would rise with absolute drift.\n")

    flux = []
    for v, _vt, r, *_ in scans:
        rs = detrend(v, smooth(r))
        m = (v > LINE_WINDOW[0]) & (v < LINE_WINDOW[1])
        flux.append(float(_trapz(rs[m], v[m])))

    flux = np.asarray(flux)
    drift = np.abs(np.asarray(drifts))
    c = float(np.corrcoef(drift, flux)[0, 1])
    clean = drift < 0.05

    print(f"   corr(|drift|, line flux)      : {c:+.2f}")
    if clean.any():
        print(f"   clean baselines (n={clean.sum():2d})        : mean flux {flux[clean].mean():+.2f}")
    if (~clean).any():
        print(f"   drifting baselines (n={(~clean).sum():2d})     : mean flux {flux[~clean].mean():+.2f}")
    verdict = ("drift is DEGRADING the data, not creating the line" if c < 0.4
               else "WARNING: line flux tracks the drift -- suspect the baseline")
    print(f"\n   VERDICT: {verdict}")


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    paths = sorted(glob.glob(os.path.join(here, "spectra", "*.npz")))
    if not paths:
        print("No spectra found in spectra/ -- nothing to check.")
        return 1

    scans = [load(p) for p in paths]
    print(f"Loaded {len(scans)} saved spectra.")

    test_sky_lock(scans)
    test_velocity_extent(scans)
    drifts = test_baseline_drift(scans)
    test_drift_vs_flux(scans, drifts)

    print("\n" + "=" * 78)
    print(" Interpretation of all four is written up in README.md.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
