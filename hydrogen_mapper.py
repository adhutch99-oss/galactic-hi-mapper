#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Galactic Plane 21-cm Hydrogen-Line Mapper
==========================================

A single-file, self-contained desktop application for a DIY radio-astronomy
telescope that maps the 21-cm neutral-hydrogen line (rest frequency
1420.405752 MHz) emitted by the galactic plane.

HARDWARE (real backend)
-----------------------
  * SDR ...... Nooelec NESDR SMArTee v2 (RTL2832U + R820T2).  The SMA port has
               an ALWAYS-ON 4.5 V bias-tee -- there is NOTHING to toggle in
               software, and the code deliberately never touches a bias-tee.
  * LNA ...... Nooelec SAWbird+ H1 (+40 dB, SAW filter @ 1420 MHz), powered by
               the SDR bias-tee.  Keep it inline at all times.
  * Antenna .. Nooelec 20 dBi parabolic mesh dish (1.4 GHz).
  * Feed ..... LMR400 (dish->LNA) + RG316 pigtail (LNA->SDR).

SOFTWARE STACK
--------------
  Python + PyQt6 (GUI) + pyrtlsdr (hardware) + numpy/scipy (DSP) +
  astropy (coordinates & velocity frames) + matplotlib (rendering).

DUAL DATA SOURCE
----------------
  Everything downstream of the acquisition layer runs identically on either
  backend.  A GUI toggle switches between:
      * RealSource -- pyrtlsdr, read_samples_async, 131072-sample chunks.
      * SimSource  -- synthetic IQ with a realistic noise floor, a modeled
                      receiver bandpass ("gain skirt"), an injected Gaussian
                      HI peak whose Doppler offset comes from a toy galactic
                      rotation model, adjustable injected SNR, and optional
                      fake RFI spikes.
  If no SDR is found on startup the app automatically falls back to Simulation
  and says so in the status bar.  Simulation is the primary development/test
  path and needs no hardware.

SELF-TEST
---------
  Run `python hydrogen_mapper.py --self-test` to exercise the entire pipeline
  head-lessly (Qt "offscreen" platform, Simulation backend):  it starts a live
  preview, logs a drift scan, records a row, and verifies that a valid line was
  written to galactic_plane_log.csv and that a point landed on both maps.
  Exit code 0 == pass, 1 == fail.

PACKAGING (PyInstaller --windowed)
----------------------------------
  astropy ships data files that a frozen build must carry.  Build with, e.g.:

      pyinstaller --windowed --name HydrogenMapper ^
          --collect-data astropy ^
          --collect-data astropy_iers_data ^
          --collect-data pyerfa ^
          --hidden-import astropy.coordinates ^
          --hidden-import scipy.signal ^
          hydrogen_mapper.py

  (Add `--hidden-import rtlsdr` only if you bundle the RTL-SDR support.)
"""

from __future__ import annotations

import csv
import math
import os
import queue
import shutil
import sys
import threading
import time
import traceback
from datetime import datetime, timedelta

import numpy as np
from scipy.signal import find_peaks, peak_widths

# --------------------------------------------------------------------------- #
#  astropy -- keep it fully offline (bundled IERS data), no network at runtime #
# --------------------------------------------------------------------------- #
from astropy.utils import iers
iers.conf.auto_download = False        # never phone home for Earth-orientation
iers.conf.auto_max_age = None          # tolerate slightly stale bundled tables

import astropy.units as u
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
from astropy.time import Time

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter("ignore", AstropyWarning)

# --------------------------------------------------------------------------- #
#  Qt                                                                          #
# --------------------------------------------------------------------------- #
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import QThread, QTimer, pyqtSignal

# --------------------------------------------------------------------------- #
#  matplotlib on the Qt canvas (no pyplot -> avoids backend surprises)         #
# --------------------------------------------------------------------------- #
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# ===========================================================================
#  CONSTANTS
# ===========================================================================
F_REST_HZ      = 1420.405752e6         # HI rest frequency (Hz)
C_KMS          = 299792.458            # speed of light (km/s)

CENTER_HZ      = 1420.00e6             # offset-tune: line sits ~0.4 MHz above DC
SAMPLE_RATE_HZ = 2.4e6                 # 2.4 MSps (fall back from 3.2 if drops)
CHUNK          = 131072                # samples per async chunk
NPERSEG        = 4096                  # Welch sub-block -> 4096 freq bins

# Peak-search window: +/-1.0 MHz around the rest line covers realistic galactic
# Doppler of roughly +/-200 km/s, and rejects band-edge / RFI artifacts outside.
SEARCH_HALF_HZ = 1.0e6

# Standard solar motion toward the LSR (Schoenrich, Binney & Dehnen 2010), km/s.
SOLAR_U, SOLAR_V, SOLAR_W = 11.1, 12.24, 7.25

# Toy flat rotation curve used ONLY by the simulator to inject a plausible line.
TOY_V0_KMS = 220.0

# Galactic rotation constants used by the target planner / analysis.
R0_KPC = 8.2                           # Sun's galactocentric distance (kpc)
V0_KMS = 220.0                         # circular speed at the Sun (km/s)

# Max HI velocity components logged per drift scan (one per arm on a sightline).
MAX_COMPONENTS = 5

# Detection smoothing: real galactic HI lines are ~10-40 km/s wide (50-200 kHz),
# far wider than one 586 Hz FFT bin.  Peak-finding runs on the spectrum smoothed
# to this width so a faint broad line beats per-bin noise (and 1-bin RFI, which
# the boxcar dilutes by the kernel length, is separately vetoed by energy
# concentration in the raw bins).
SMOOTH_DETECT_HZ = 50e3

# Blocks averaged for the cold-sky baseline.  The OFF spectrum divides every
# later scan, so its noise floor imprints on all of them: a single display
# frame (~0.15 s of data) is nowhere near enough.  150 blocks ~ 8 s of raw
# data on hardware, ~7x less baseline noise than the old EMA snapshot.
BASELINE_BLOCKS = 150

# Sidereal hours -> solar hours conversion (a sidereal day is shorter).
SIDEREAL_TO_SOLAR = 23.9344696 / 24.0

# Valid R820T2 tuner gains (dB).  Default is high-but-not-max: with +40 dB of
# LNA ahead of it, slamming the SDR to 49.6 dB saturates and adds noise.
RTL_GAINS_DB = [0.0, 0.9, 1.4, 2.7, 3.7, 7.7, 8.7, 12.5, 14.4, 15.7, 16.6,
                19.7, 20.7, 22.9, 25.4, 28.0, 29.7, 32.8, 33.8, 36.4, 37.2,
                38.6, 40.2, 42.1, 43.4, 43.9, 44.5, 48.0, 49.6]
DEFAULT_GAIN_DB = 40.2

# Map geometry.
L_BINS = 360                           # galactic longitude 0..360 deg, 1 deg/bin
B_MIN, B_MAX = -20, 20                 # galactic latitude range (deg)
B_BINS = B_MAX - B_MIN                 # 1 deg/bin

# The original spec columns, plus: Source (Sim/Hardware) so simulated and real
# rows are never confused, and Component_Index / N_Components so several HI
# velocity components (multiple arms on one sightline) can share a scan.
# SNR_sigma / Line_FWHM_kHz / Integration_s / Spectrum_File were added so the
# analysis can filter junk fallback rows and rebuild full line profiles.
CSV_HEADER = [
    "Timestamp_UTC", "Source", "Lat_GPS", "Lon_GPS", "Dish_Azimuth",
    "Dish_Altitude", "Galactic_Longitude_L", "Galactic_Latitude_B",
    "Peak_Frequency_MHz", "Peak_Power_dB", "Doppler_Velocity_Topo_km_s",
    "Doppler_Velocity_LSR_km_s", "Component_Index", "N_Components",
    "SNR_sigma", "Line_FWHM_kHz", "Integration_s", "Spectrum_File",
]

# Log file lives next to this script (so it survives a PyInstaller onefile run
# via the executable's directory).
def _app_dir() -> str:
    if getattr(sys, "frozen", False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

CSV_PATH = os.path.join(_app_dir(), "galactic_plane_log.csv")

# Every recorded scan also saves its full averaged spectrum here (.npz), so the
# l-v diagram can be rebuilt as a true 2-D intensity image, not just peaks.
SPECTRA_DIR = os.path.join(_app_dir(), "spectra")

# Longitude ranges that never transit due south of the fixed Long Island site
# (l ~ 90-165 crosses north of zenith; l ~ 275-335 stays below the horizon).
UNREACHABLE_L = ((90, 165), (275, 335))


# ===========================================================================
#  DSP HELPERS
# ===========================================================================
def compute_psd(iq: np.ndarray, fs: float, center: float,
                nperseg: int = NPERSEG):
    """Welch-style averaged two-sided PSD of a complex IQ chunk.

    Returns (freqs_hz, power_db) where freqs_hz is the absolute RF axis
    (center +/- fs/2) and power_db is 10*log10 of the averaged periodogram.
    """
    n = len(iq)
    if n < nperseg:
        nperseg = 1 << int(math.log2(max(n, 2)))
    win = np.hanning(nperseg)
    win_power = np.sum(win ** 2)
    step = nperseg // 2
    acc = np.zeros(nperseg, dtype=np.float64)
    count = 0
    for start in range(0, n - nperseg + 1, step):
        seg = iq[start:start + nperseg] * win
        sp = np.fft.fftshift(np.fft.fft(seg))
        acc += np.abs(sp) ** 2
        count += 1
    if count == 0:                     # pathological tiny chunk
        seg = iq[:nperseg] * win[:len(iq[:nperseg])]
        acc = np.abs(np.fft.fftshift(np.fft.fft(seg, n=nperseg))) ** 2
        count = 1
    psd = acc / count
    psd /= (fs * win_power)            # power spectral density normalisation
    power_db = 10.0 * np.log10(psd + 1e-20)
    freqs = center + np.fft.fftshift(np.fft.fftfreq(nperseg, d=1.0 / fs))
    return freqs.astype(np.float64), power_db.astype(np.float64)


def notch_dc(freqs: np.ndarray, power_db: np.ndarray, center: float,
             half_bins: int = 3) -> np.ndarray:
    """Interpolate across the RTL-SDR DC spike (baseband 0 Hz == `center`).

    The line sits ~0.4 MHz off DC, but the +/-1 MHz search window still spans
    DC, so the residual spike must be removed before peak-finding.
    """
    out = power_db.copy()
    k = int(np.argmin(np.abs(freqs - center)))
    lo = max(k - half_bins, 0)
    hi = min(k + half_bins, len(out) - 1)
    left = out[lo - 1] if lo - 1 >= 0 else out[hi + 1]
    right = out[hi + 1] if hi + 1 < len(out) else out[lo - 1]
    out[lo:hi + 1] = np.linspace(left, right, hi + 1 - lo)
    return out


# ===========================================================================
#  ASTROPHYSICS
# ===========================================================================
def astro_solve(lat_deg: float, lon_deg: float, alt_deg: float,
                az_deg: float = 180.0, when: Time | None = None) -> dict:
    """Full coordinate + velocity-frame solution for the current pointing.

    Returns a dict with galactic (l, b), UTC/JD/LST, the SkyCoord, and the
    correction (km/s) that converts a topocentric radial velocity to the LSR:

        v_lsr = v_topo + lsr_corr

    where lsr_corr = barycentric_correction + solar_motion_projection.
    """
    t = when if when is not None else Time.now()
    loc = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg, height=100 * u.m)

    altaz_frame = AltAz(obstime=t, location=loc)
    coo = SkyCoord(az=az_deg * u.deg, alt=alt_deg * u.deg, frame=altaz_frame)
    icrs = coo.icrs
    gal = coo.galactic
    l_deg = float(gal.l.deg)
    b_deg = float(gal.b.deg)

    # Barycentric (Earth motion) correction -- add to observed radial velocity.
    # Use a bare ra/dec SkyCoord so it carries no obstime/location frame
    # attributes (those would clash with the explicit arguments below).
    sc_plain = SkyCoord(ra=icrs.ra, dec=icrs.dec, frame="icrs")
    bary_kms = float(sc_plain.radial_velocity_correction(
        kind="barycentric", obstime=t, location=loc).to(u.km / u.s).value)

    # Projection of the Sun's peculiar motion onto the line of sight.
    lr, br = math.radians(l_deg), math.radians(b_deg)
    solar_kms = (SOLAR_U * math.cos(lr) * math.cos(br)
                 + SOLAR_V * math.sin(lr) * math.cos(br)
                 + SOLAR_W * math.sin(br))

    lsr_corr = bary_kms + solar_kms

    try:
        lst = t.sidereal_time("apparent", longitude=lon_deg * u.deg).to_string(unit=u.hour)
    except Exception:
        lst = "n/a"

    return {
        "time": t, "utc": t.isot, "jd": float(t.jd), "lst": lst,
        "l": l_deg, "b": b_deg, "skycoord": icrs,
        "bary_kms": bary_kms, "solar_kms": solar_kms, "lsr_corr": lsr_corr,
    }


def toy_injected_v_lsr(l_deg: float, b_deg: float) -> float:
    """Toy flat-rotation-curve LSR velocity the simulator injects for (l, b)."""
    return TOY_V0_KMS * math.sin(math.radians(l_deg)) * math.cos(math.radians(b_deg))


def plan_pointing(l_deg: float, lat_deg: float, lon_deg: float) -> dict:
    """Where/when to point the fixed south beam to catch galactic longitude l.

    Returns the meridian-transit altitude (dish faces due south at transit),
    the source declination, whether it is reachable due south, and the local
    clock time of the next southern transit.
    """
    gc = SkyCoord(l=l_deg * u.deg, b=0 * u.deg, frame="galactic")
    icrs = gc.icrs
    ra_h = float(icrs.ra.hour)
    dec = float(icrs.dec.deg)

    # At upper transit an object due south culminates at alt = 90 - |lat - dec|.
    alt_transit = 90.0 - abs(lat_deg - dec)
    # It transits to the SOUTH (az=180) only if its declination is below the
    # observer's latitude; otherwise it crosses north of the zenith.
    transits_south = dec <= lat_deg
    observable = transits_south and alt_transit > 0.0

    # Local clock time of the next transit: LST must equal the source RA.
    now = Time.now()
    lst_now = float(now.sidereal_time("apparent", longitude=lon_deg * u.deg).hour)
    dh_sidereal = (ra_h - lst_now) % 24.0
    dh_solar = dh_sidereal * SIDEREAL_TO_SOLAR
    t_local = datetime.now() + timedelta(hours=dh_solar)

    return {
        "ra_h": ra_h, "dec": dec, "alt_transit": alt_transit,
        "transits_south": transits_south, "observable": observable,
        "hours_to_transit": dh_solar, "transit_local": t_local,
    }


# ===========================================================================
#  DATA SOURCES  (one interface, two interchangeable backends)
# ===========================================================================
class DataSource:
    """Common interface: open() -> repeated read_iq() -> close()."""

    name = "base"

    def __init__(self, sample_rate=SAMPLE_RATE_HZ, center=CENTER_HZ,
                 gain_db=DEFAULT_GAIN_DB, chunk=CHUNK):
        self.sample_rate = sample_rate
        self.center = center
        self.gain_db = gain_db
        self.chunk = chunk

    def open(self):     raise NotImplementedError
    def read_iq(self) -> np.ndarray: raise NotImplementedError
    def close(self):    pass


class SimSource(DataSource):
    """Synthetic IQ generator with a physically-motivated PSD.

    The line is injected as *narrowband noise* (not a coherent tone) because
    real HI emission is broadband thermal radiation.  Its baseband frequency is
    driven by `f_peak_baseband`, which the GUI updates from the toy rotation
    model so that the analysis path recovers the injected v_lsr round-trip.
    """

    name = "Simulation"

    def __init__(self, **kw):
        super().__init__(**kw)
        self._rng = np.random.default_rng()
        # Default: line at rest (no Doppler) -> baseband offset of the rest line.
        self.f_peak_baseband = F_REST_HZ - self.center
        self.snr_db = 12.0             # injected peak-over-floor (dB)
        self.line_sigma_hz = 80e3      # ~17 km/s velocity width
        self.rfi_enabled = False
        # Precompute the baseband frequency axis for a full chunk.
        self._f = np.fft.fftshift(np.fft.fftfreq(self.chunk, d=1.0 / self.sample_rate))
        self._bp = self._bandpass(self._f)

    # ---- receiver bandpass / "gain skirt" -------------------------------- #
    def _bandpass(self, f):
        """Linear-amplitude passband shape (SAW filter + tuner roll-off)."""
        # Broad SAW passband hump plus gentle skirts toward the band edges.
        passband = 0.35 + 0.65 * np.exp(-0.5 * (f / 0.85e6) ** 2)
        tilt = 1.0 + 0.05 * (f / (self.sample_rate / 2))   # slight asymmetry
        return passband * tilt

    def open(self):
        # nothing to power up in simulation
        return True

    def read_iq(self) -> np.ndarray:
        n = self.chunk
        rng = self._rng
        f = self._f
        bp = self._bp

        # Gain scales the floor a touch (visual realism only).
        gain_lin = 10 ** ((self.gain_db - DEFAULT_GAIN_DB) / 40.0)

        # 1) White complex thermal noise -> frequency domain.
        base = (rng.standard_normal(n) + 1j * rng.standard_normal(n))
        S = np.fft.fftshift(np.fft.fft(base))

        # 2) Colour it with the receiver bandpass.
        S *= bp * gain_lin

        # 3) Inject the HI line as narrowband noise (Gaussian power profile).
        line_shape = np.exp(-0.5 * ((f - self.f_peak_baseband) / self.line_sigma_hz) ** 2)
        line_gain = np.sqrt(n) * bp * gain_lin * math.sqrt(10 ** (self.snr_db / 10.0))
        line_noise = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / math.sqrt(2)
        S += line_gain * line_shape * line_noise

        # 4) Modest residual DC spike (single-bin RTL artifact) for realism.
        dc = int(np.argmin(np.abs(f)))
        S[dc] += 8.0 * math.sqrt(n) * float(np.max(bp))

        # 5) Optional fake RFI: sharp coherent spikes at random baseband freqs.
        if self.rfi_enabled:
            for _ in range(2):
                rf = rng.uniform(-0.45 * self.sample_rate, 0.45 * self.sample_rate)
                idx = int(np.argmin(np.abs(f - rf)))
                S[idx] += 40.0 * math.sqrt(n) * float(np.max(bp))

        iq = np.fft.ifft(np.fft.ifftshift(S))
        return iq.astype(np.complex64)


class RealSource(DataSource):
    """pyrtlsdr backend using read_samples_async into a bounded queue.

    The async callback (worker of librtlsdr) pushes chunks; read_iq() pulls the
    newest one so the DSP thread sees the same pull-based API as SimSource.
    NOTE: the SMArTee v2 bias-tee is ALWAYS ON in hardware -- we never toggle it.
    """

    name = "Hardware"

    def __init__(self, **kw):
        super().__init__(**kw)
        self.sdr = None
        self._q: queue.Queue = queue.Queue(maxsize=8)
        self._async_thread = None
        self._closing = False

    def open(self):
        from rtlsdr import RtlSdr           # imported lazily: only if hardware used
        self.sdr = RtlSdr()
        self.sdr.sample_rate = self.sample_rate
        self.sdr.center_freq = self.center  # offset tune @ 1420.00 MHz
        try:
            self.sdr.gain = self.gain_db
        except Exception:
            self.sdr.gain = "auto"
        self._closing = False
        self._async_thread = threading.Thread(target=self._async_loop, daemon=True)
        self._async_thread.start()
        return True

    def _async_loop(self):
        try:
            self.sdr.read_samples_async(self._on_samples, self.chunk)
        except Exception:
            pass                            # cancelled on close

    def _on_samples(self, samples, _ctx):
        if self._closing:
            return
        try:
            self._q.put_nowait(np.asarray(samples, dtype=np.complex64))
        except queue.Full:
            try:                            # drop oldest, keep freshest
                self._q.get_nowait()
                self._q.put_nowait(np.asarray(samples, dtype=np.complex64))
            except queue.Empty:
                pass

    def set_gain(self, gain_db):
        self.gain_db = gain_db
        if self.sdr is not None:
            try:
                self.sdr.gain = gain_db
            except Exception:
                pass

    def read_iq(self) -> np.ndarray:
        return self._q.get(timeout=5.0)

    def close(self):
        self._closing = True
        if self.sdr is not None:
            try:
                self.sdr.cancel_read_async()
            except Exception:
                pass
        if self._async_thread is not None:
            self._async_thread.join(timeout=2.0)
        if self.sdr is not None:
            try:
                self.sdr.close()
            except Exception:
                pass
        self.sdr = None


def make_source(kind: str, gain_db: float) -> DataSource:
    if kind == "Hardware":
        return RealSource(gain_db=gain_db)
    return SimSource(gain_db=gain_db)


def hardware_available() -> bool:
    """True only if an RTL-SDR can actually be opened right now."""
    try:
        from rtlsdr import RtlSdr
        sdr = RtlSdr()
        sdr.close()
        return True
    except Exception:
        return False


# ===========================================================================
#  ACQUISITION WORKER  (runs the read/DSP loop off the GUI thread)
# ===========================================================================
class AcquisitionWorker(QThread):
    """Owns the DataSource, produces notched PSDs, and marshals them to the GUI
    via a queued signal.  It NEVER touches Qt widgets or matplotlib directly."""

    spectrum_ready = pyqtSignal(object, object)   # (freqs_hz, power_db)
    status = pyqtSignal(str)
    failed = pyqtSignal(str)

    def __init__(self, source: DataSource, parent=None):
        super().__init__(parent)
        self.source = source
        self._running = False

    def run(self):
        self._running = True
        try:
            self.source.open()
        except Exception as exc:
            self.failed.emit(f"Source open failed: {exc}")
            return

        self.status.emit(f"Streaming from {self.source.name}")
        is_sim = isinstance(self.source, SimSource)
        while self._running:
            try:
                iq = self.source.read_iq()
            except Exception as exc:
                if self._running:
                    self.failed.emit(f"Read error: {exc}")
                break
            freqs, power_db = compute_psd(iq, self.source.sample_rate,
                                          self.source.center)
            power_db = notch_dc(freqs, power_db, self.source.center)
            self.spectrum_ready.emit(freqs, power_db)
            if is_sim:
                self.msleep(40)            # ~25 fps for a smooth simulated view
        try:
            self.source.close()
        except Exception:
            pass

    def stop(self):
        self._running = False
        self.wait(4000)


# ===========================================================================
#  THEME (dark, night-observing friendly) + MATPLOTLIB CANVAS WRAPPERS
# ===========================================================================
PLOT_BG   = "#15181c"      # figure background
PLOT_AX   = "#0d0f12"      # axes background
PLOT_FG   = "#d0d7e0"      # text / ticks
PLOT_GRID = "#3a3f46"      # spines / outlines


def apply_dark_theme(app: QtWidgets.QApplication):
    """App-wide dark Fusion palette + a readable default font."""
    app.setStyle("Fusion")
    app.setFont(QtGui.QFont("Segoe UI", 10))
    C = QtGui.QColor
    R = QtGui.QPalette.ColorRole
    G = QtGui.QPalette.ColorGroup
    p = QtGui.QPalette()
    p.setColor(R.Window, C("#1e2228"))
    p.setColor(R.WindowText, C("#d8dee6"))
    p.setColor(R.Base, C("#15181c"))
    p.setColor(R.AlternateBase, C("#1b1f24"))
    p.setColor(R.Text, C("#d8dee6"))
    p.setColor(R.Button, C("#2a2f36"))
    p.setColor(R.ButtonText, C("#d8dee6"))
    p.setColor(R.Highlight, C("#2f6db3"))
    p.setColor(R.HighlightedText, C("#ffffff"))
    p.setColor(R.ToolTipBase, C("#2a2f36"))
    p.setColor(R.ToolTipText, C("#d8dee6"))
    p.setColor(R.PlaceholderText, C("#7a828c"))
    p.setColor(R.Link, C("#4da3ff"))
    for role in (R.Text, R.ButtonText, R.WindowText):
        p.setColor(G.Disabled, role, C("#6a7078"))
    app.setPalette(p)


def style_axis(ax):
    """Match a matplotlib axes to the dark GUI theme."""
    ax.set_facecolor(PLOT_AX)
    for sp in ax.spines.values():
        sp.set_color(PLOT_GRID)
    ax.tick_params(colors=PLOT_FG, labelsize=8)
    ax.xaxis.label.set_color(PLOT_FG)
    ax.yaxis.label.set_color(PLOT_FG)
    ax.title.set_color(PLOT_FG)


def style_colorbar(cbar):
    cbar.ax.tick_params(colors=PLOT_FG, labelsize=7)
    cbar.outline.set_edgecolor(PLOT_GRID)


def freq_mhz_to_vlsr_axis(f_mhz):
    """Topocentric radio-convention velocity (km/s) for the secondary axis."""
    return C_KMS * (F_REST_HZ - np.asarray(f_mhz) * 1e6) / F_REST_HZ


def vlsr_axis_to_freq_mhz(v_kms):
    return F_REST_HZ * (1.0 - np.asarray(v_kms) / C_KMS) / 1e6


class MplCanvas(FigureCanvas):
    def __init__(self, height=2.6):
        self.fig = Figure(figsize=(5, height), tight_layout=True,
                          facecolor=PLOT_BG)
        super().__init__(self.fig)
        self.ax = self.fig.add_subplot(111)
        style_axis(self.ax)


# ===========================================================================
#  MAIN WINDOW
# ===========================================================================
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, force_source: str | None = None):
        super().__init__()
        self.setWindowTitle("Galactic Plane 21-cm HI Mapper")
        self.resize(1280, 820)

        # ---- runtime state ------------------------------------------------ #
        self.worker: AcquisitionWorker | None = None
        self.source: DataSource | None = None

        self.freqs = None                 # RF axis (Hz) of the current spectrum
        self.last_spectrum = None         # most recent raw PSD (dB, diagnostics)
        self.ema_lin = None               # smoothed LINEAR PSD (display/baseline)
        self.baseline_lin = None          # captured cold-sky passband (LINEAR)
        self._bl_sum = None               # baseline accumulator (linear power)
        self._bl_count = 0
        self._bl_target = 0               # blocks still wanted; 0 = idle

        self.is_logging = False
        self.log_sum = None               # running LINEAR-power sum (radiometry:
        self.log_count = 0                # average power, never average dB)
        self._log_t0: float | None = None # wall-clock start of the integration
        self._log_last_data_t: float | None = None  # wall-clock of last block
        self._force_rescale = False       # one-shot y rescale past a frozen axis

        self.spectrum_count = 0           # diagnostics / self-test
        self.integration_count = 0

        self.records: list[dict] = []     # every logged row
        self.lb_map = np.full((B_BINS, L_BINS), np.nan)  # [b, l] -> peak dB
        self.lv_l: list[float] = []       # longitude-velocity scatter data
        self.lv_v: list[float] = []
        self.lv_p: list[float] = []

        # injection parameters the simulator reads (thread-safe: plain floats)
        self._inj = {"l": 0.0, "b": 0.0, "v_lsr": 0.0, "lsr_corr": 0.0}

        # Optional fixed observation time (used only by the automated self-test
        # so the sky geometry is deterministic); None == use the real clock.
        self._fixed_time: Time | None = None
        self._planned_alt: float | None = None   # last planner result

        self._forced_source = force_source

        self._build_ui()
        self._decide_default_source()
        self._ensure_csv_schema()
        self._load_existing_log()
        self._refresh_maps()

        # Periodically refresh the simulated injection so it tracks time/pointing.
        self._inj_timer = QTimer(self)
        self._inj_timer.timeout.connect(self.update_injection)
        self._inj_timer.start(2000)

    # ------------------------------------------------------------------ UI --
    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        root = QtWidgets.QHBoxLayout(central)

        # ---------------- LEFT: control panel (scrollable) ---------------- #
        panel = QtWidgets.QWidget()
        pl = QtWidgets.QVBoxLayout(panel)
        pl.setSpacing(6)

        # Big ON-PLANE banner: the one judgment call the whole workflow hinges
        # on (b ~ 0 = transit = record now) gets a traffic-light indicator.
        self.plane_banner = QtWidgets.QLabel("Pointing: computing…")
        self.plane_banner.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.plane_banner.setWordWrap(True)
        self.plane_banner.setMinimumHeight(48)
        self._set_banner("#4a4f57", "Pointing: computing…")
        pl.addWidget(self.plane_banner)

        # GPS
        gps_box = QtWidgets.QGroupBox("GPS Location")
        gps_form = QtWidgets.QFormLayout(gps_box)
        self.lat_spin = QtWidgets.QDoubleSpinBox()
        self.lat_spin.setRange(-90.0, 90.0); self.lat_spin.setDecimals(5)
        self.lat_spin.setValue(40.70); self.lat_spin.setSuffix(" °")   # default site: Long Island, NY (rounded)
        self.lon_spin = QtWidgets.QDoubleSpinBox()
        self.lon_spin.setRange(-180.0, 180.0); self.lon_spin.setDecimals(5)
        self.lon_spin.setValue(-73.60); self.lon_spin.setSuffix(" °")   # W = negative (rounded)
        gps_form.addRow("Latitude:", self.lat_spin)
        gps_form.addRow("Longitude:", self.lon_spin)
        pl.addWidget(gps_box)

        # Alignment
        al_box = QtWidgets.QGroupBox("Dish Alignment")
        al_form = QtWidgets.QFormLayout(al_box)
        az_label = QtWidgets.QLabel("180  (South)")
        az_label.setStyleSheet("font-weight: bold;")
        self.alt_spin = QtWidgets.QDoubleSpinBox()
        self.alt_spin.setRange(0.0, 90.0); self.alt_spin.setDecimals(2)
        self.alt_spin.setValue(45.0); self.alt_spin.setSuffix(" °")
        al_form.addRow("Azimuth:", az_label)
        al_form.addRow("Altitude (deg):", self.alt_spin)
        pl.addWidget(al_box)
        self.alt_spin.valueChanged.connect(lambda *_: self.update_injection())
        self.lat_spin.valueChanged.connect(lambda *_: self.update_injection())
        self.lon_spin.valueChanged.connect(lambda *_: self.update_injection())

        # Source + gain
        src_box = QtWidgets.QGroupBox("Data Source")
        src_form = QtWidgets.QFormLayout(src_box)
        self.source_combo = QtWidgets.QComboBox()
        self.source_combo.addItems(["Simulation", "Hardware"])
        self.gain_combo = QtWidgets.QComboBox()
        for g in RTL_GAINS_DB:
            self.gain_combo.addItem(f"{g:.1f} dB", g)
        self.gain_combo.setCurrentIndex(RTL_GAINS_DB.index(DEFAULT_GAIN_DB))
        self.gain_combo.currentIndexChanged.connect(self._on_gain_changed)
        src_form.addRow("Backend:", self.source_combo)
        src_form.addRow("SDR Gain:", self.gain_combo)
        pl.addWidget(src_box)

        # Simulation controls
        sim_box = QtWidgets.QGroupBox("Simulation Injection")
        sim_form = QtWidgets.QFormLayout(sim_box)
        self.snr_spin = QtWidgets.QDoubleSpinBox()
        self.snr_spin.setRange(0.0, 40.0); self.snr_spin.setValue(12.0)
        self.snr_spin.setSuffix(" dB")
        self.snr_spin.valueChanged.connect(self._on_sim_changed)
        self.rfi_check = QtWidgets.QCheckBox("Inject fake RFI spikes")
        self.rfi_check.stateChanged.connect(self._on_sim_changed)
        sim_form.addRow("Injected SNR:", self.snr_spin)
        sim_form.addRow(self.rfi_check)
        pl.addWidget(sim_box)

        # Target planner: "I want longitude l -> which altitude & when?"
        plan_box = QtWidgets.QGroupBox("Target Planner (l -> altitude & time)")
        plan_v = QtWidgets.QVBoxLayout(plan_box)
        plan_form = QtWidgets.QFormLayout()
        self.target_l_spin = QtWidgets.QDoubleSpinBox()
        self.target_l_spin.setRange(0.0, 360.0); self.target_l_spin.setDecimals(1)
        self.target_l_spin.setValue(40.0); self.target_l_spin.setSuffix(" °")
        plan_form.addRow("Galactic l:", self.target_l_spin)
        plan_v.addLayout(plan_form)
        plan_btns = QtWidgets.QHBoxLayout()
        self.btn_plan = QtWidgets.QPushButton("Find pointing")
        self.btn_plan_apply = QtWidgets.QPushButton("Apply altitude")
        plan_btns.addWidget(self.btn_plan)
        plan_btns.addWidget(self.btn_plan_apply)
        plan_v.addLayout(plan_btns)
        self.btn_tonight = QtWidgets.QPushButton("What can I see tonight?")
        plan_v.addWidget(self.btn_tonight)
        self.plan_label = QtWidgets.QLabel("Enter a longitude and click "
                                           "“Find pointing”.")
        self.plan_label.setWordWrap(True)
        plan_v.addWidget(self.plan_label)
        pl.addWidget(plan_box)
        self.btn_plan.clicked.connect(self.plan_target)
        self.btn_plan_apply.clicked.connect(self.apply_planned_altitude)
        self.btn_tonight.clicked.connect(self.whats_up_tonight)

        # Actions
        act_box = QtWidgets.QGroupBox("Actions")
        av = QtWidgets.QVBoxLayout(act_box)
        self.btn_preview = QtWidgets.QPushButton("START LIVE PREVIEW")
        self.btn_log = QtWidgets.QPushButton("LOG DRIFT SCAN")
        self.btn_record = QtWidgets.QPushButton("STOP && RECORD ROW")
        self.btn_baseline = QtWidgets.QPushButton("Capture Cold Sky Baseline")
        self.chk_subtract = QtWidgets.QCheckBox("Subtract baseline from spectrum")
        for b in (self.btn_preview, self.btn_log, self.btn_record, self.btn_baseline):
            b.setMinimumHeight(34)
            av.addWidget(b)
        av.addWidget(self.chk_subtract)
        self.chk_freeze_y = QtWidgets.QCheckBox("Freeze spectrum y-axis")
        av.addWidget(self.chk_freeze_y)
        self.integ_label = QtWidgets.QLabel("Not integrating.")
        self.integ_label.setStyleSheet("color: #9aa3ad;")
        av.addWidget(self.integ_label)
        note = QtWidgets.QLabel("Tip: capture the baseline while pointed OFF the "
                                "galactic plane (cold sky).")
        note.setWordWrap(True)
        note.setStyleSheet("color: #9aa3ad;")
        av.addWidget(note)
        pl.addWidget(act_box)

        self.btn_preview.clicked.connect(self.start_preview)
        self.btn_log.clicked.connect(self.start_logging)
        self.btn_record.clicked.connect(self.stop_and_record)
        self.btn_baseline.clicked.connect(self.capture_baseline)
        self.chk_subtract.stateChanged.connect(self._on_subtract_toggled)
        self.chk_freeze_y.stateChanged.connect(lambda *_: self._redraw_spectrum())

        self.btn_log.setEnabled(False)
        self.btn_record.setEnabled(False)
        self.btn_baseline.setEnabled(False)

        pl.addStretch(1)

        self.info_label = QtWidgets.QLabel("")
        self.info_label.setWordWrap(True)
        # Monospace: these are the numbers you read from across the yard.
        self.info_label.setStyleSheet(
            "font-family: Consolas, monospace; font-size: 12px;")
        pl.addWidget(self.info_label)

        # Put the (tall) control panel in a fixed-width scroll area so nothing
        # is clipped on shorter screens.
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(panel)
        scroll.setWidgetResizable(True)
        scroll.setFixedWidth(360)
        scroll.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        root.addWidget(scroll)

        # ---------------- RIGHT: three stacked canvases ---------------- #
        right = QtWidgets.QWidget()
        rl = QtWidgets.QVBoxLayout(right)
        rl.setSpacing(4)

        self.c_spec = MplCanvas(height=2.4)
        self.c_map = MplCanvas(height=2.4)
        self.c_lv = MplCanvas(height=2.4)
        rl.addWidget(self.c_spec)
        rl.addWidget(self.c_map)
        rl.addWidget(self.c_lv)
        root.addWidget(right, 1)

        self._init_plots()

        # status bar
        self.status_lbl = QtWidgets.QLabel("Ready.")
        self.statusBar().addWidget(self.status_lbl)

    def _init_plots(self):
        # --- 1D spectrum (velocity axis on top: this IS a Doppler instrument) ---
        ax = self.c_spec.ax
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel("Power (dB)")
        (self.spec_line,) = ax.plot([], [], lw=0.6, color="#4da3ff", alpha=0.55)
        # Bold overlay: the spectrum smoothed to HI line width.  A real (broad)
        # hydrogen line stands out here long before it is visible per-bin.
        (self.spec_line_sm,) = ax.plot([], [], lw=1.8, color="#ffd166",
                                       label="smoothed (50 kHz)")
        lo = (F_REST_HZ - SEARCH_HALF_HZ) / 1e6
        hi = (F_REST_HZ + SEARCH_HALF_HZ) / 1e6
        self.search_span = ax.axvspan(lo, hi, color="orange", alpha=0.12,
                                      label="search window")
        self.rest_line = ax.axvline(F_REST_HZ / 1e6, color="#4caf50", ls="--",
                                    lw=0.8, label="rest 1420.4058 (v=0)")
        (self.peak_marker,) = ax.plot([], [], "rv", ms=9, label="peak")
        ax.legend(loc="upper right", fontsize=7, facecolor=PLOT_AX,
                  edgecolor=PLOT_GRID, labelcolor=PLOT_FG)
        ax.set_xlim((CENTER_HZ - SAMPLE_RATE_HZ / 2) / 1e6,
                    (CENTER_HZ + SAMPLE_RATE_HZ / 2) / 1e6)
        self.spec_vax = ax.secondary_xaxis(
            "top", functions=(freq_mhz_to_vlsr_axis, vlsr_axis_to_freq_mhz))
        self.spec_vax.set_xlabel("velocity (km/s, radio conv.; + = receding)",
                                 fontsize=8, color=PLOT_FG)
        self.spec_vax.tick_params(colors=PLOT_FG, labelsize=7)
        for sp in self.spec_vax.spines.values():
            sp.set_color(PLOT_GRID)
        self.c_spec.draw()

        # --- 2D l-b heatmap ---
        axm = self.c_map.ax
        axm.set_title("Galactic Plane Map (Peak Power)", fontsize=9)
        axm.set_xlabel("Galactic Longitude l (deg)")
        axm.set_ylabel("Galactic Latitude b (deg)")
        self.map_im = axm.imshow(self.lb_map, origin="lower", aspect="auto",
                                 extent=[0, 360, B_MIN, B_MAX], cmap="viridis",
                                 vmin=-90, vmax=-20, interpolation="nearest")
        self.map_im.cmap.set_bad(color="#101010")
        for lo_l, hi_l in UNREACHABLE_L:       # sky this site can never reach
            axm.axvspan(lo_l, hi_l, color="#666666", alpha=0.18, lw=0)
        self.map_cbar = self.c_map.fig.colorbar(self.map_im, ax=axm,
                                                fraction=0.046, pad=0.02)
        self.map_cbar.set_label("Peak Power (dB)", fontsize=8, color=PLOT_FG)
        style_colorbar(self.map_cbar)
        self.c_map.draw()

        # --- longitude-velocity diagram (the money plot) ---
        axv = self.c_lv.ax
        axv.set_title("Longitude–Velocity Diagram (Rotation Curve)", fontsize=9)
        axv.set_xlabel("Galactic Longitude l (deg)")
        axv.set_ylabel("v_LSR (km/s)")
        axv.set_xlim(0, 360)
        axv.set_ylim(-250, 250)
        axv.axhline(0, color="gray", lw=0.5)
        for lo_l, hi_l in UNREACHABLE_L:
            axv.axvspan(lo_l, hi_l, color="#666666", alpha=0.18, lw=0)
        self.lv_scatter = axv.scatter([], [], c=[], cmap="plasma",
                                      vmin=-90, vmax=-20, s=40,
                                      edgecolors="#cccccc", linewidths=0.3)
        self.lv_cbar = self.c_lv.fig.colorbar(self.lv_scatter, ax=axv,
                                              fraction=0.046, pad=0.02)
        self.lv_cbar.set_label("Peak Power (dB)", fontsize=8, color=PLOT_FG)
        style_colorbar(self.lv_cbar)
        self.c_lv.draw()

    # ---------------------------------------------------- source selection --
    def _decide_default_source(self):
        if self._forced_source == "sim":
            self.source_combo.setCurrentText("Simulation")
            self._set_status("Simulation forced (self-test / no hardware).")
            return
        if hardware_available():
            self.source_combo.setCurrentText("Hardware")
            self._set_status("RTL-SDR detected. Hardware backend available.")
        else:
            self.source_combo.setCurrentText("Simulation")
            self.source_combo.setEnabled(True)
            self._set_status("No SDR found — falling back to Simulation mode.")

    # ------------------------------------------------------ persistence -----
    def _ensure_csv_schema(self):
        """Migrate an older-format log (pre Source/Component columns) in place,
        keeping a .bak backup, so appends never produce a ragged CSV."""
        if not os.path.exists(CSV_PATH):
            return
        try:
            with open(CSV_PATH, newline="") as fh:
                reader = csv.DictReader(fh)
                existing = reader.fieldnames or []
                old_rows = list(reader)
            if existing == CSV_HEADER:
                return
            shutil.copy2(CSV_PATH, CSV_PATH + ".bak")
            with open(CSV_PATH, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=CSV_HEADER)
                writer.writeheader()
                for r in old_rows:
                    new = {k: r.get(k, "") for k in CSV_HEADER}
                    if not new.get("Source"):
                        new["Source"] = "Unknown"
                    if not new.get("Component_Index"):
                        new["Component_Index"] = "1"
                    if not new.get("N_Components"):
                        new["N_Components"] = "1"
                    writer.writerow(new)
            self._set_status("Migrated log to new schema "
                             f"(backup: {os.path.basename(CSV_PATH)}.bak).")
        except Exception as exc:
            self._set_status(f"Log schema migration skipped: {exc}")

    def _load_existing_log(self):
        if not os.path.exists(CSV_PATH):
            return
        try:
            with open(CSV_PATH, newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    try:
                        rec = {
                            "l": float(row["Galactic_Longitude_L"]),
                            "b": float(row["Galactic_Latitude_B"]),
                            "power": float(row["Peak_Power_dB"]),
                            "v_lsr": float(row["Doppler_Velocity_LSR_km_s"]),
                        }
                    except (KeyError, ValueError):
                        continue
                    self.records.append(rec)
                    self._add_point_to_maps(rec["l"], rec["b"], rec["power"],
                                            rec["v_lsr"])
            self._set_status(f"Loaded {len(self.records)} prior record(s) "
                             f"from {os.path.basename(CSV_PATH)}.")
        except Exception as exc:
            self._set_status(f"Could not read existing log: {exc}")

    # ------------------------------------------------------ gain / sim -----
    def _on_gain_changed(self, *_):
        g = self.gain_combo.currentData()
        if isinstance(self.source, RealSource):
            self.source.set_gain(g)
        elif self.source is not None:
            self.source.gain_db = g

    def _on_sim_changed(self, *_):
        if isinstance(self.source, SimSource):
            self.source.snr_db = self.snr_spin.value()
            self.source.rfi_enabled = self.rfi_check.isChecked()

    # ------------------------------------------------- injection tracking --
    def _set_banner(self, color: str, text: str):
        self.plane_banner.setText(text)
        self.plane_banner.setStyleSheet(
            f"background: {color}; color: white; font-size: 13px; "
            f"font-weight: bold; border-radius: 6px; padding: 6px;")

    def update_injection(self):
        """Recompute the toy-model injection and push it to a live SimSource."""
        try:
            sol = astro_solve(self.lat_spin.value(), self.lon_spin.value(),
                              self.alt_spin.value(), when=self._fixed_time)
        except Exception as exc:
            self._set_status(f"Astro solve failed: {exc}")
            return
        v_lsr = toy_injected_v_lsr(sol["l"], sol["b"])
        v_topo = v_lsr - sol["lsr_corr"]
        f_obs = F_REST_HZ * (1.0 - v_topo / C_KMS)
        self._inj = {"l": sol["l"], "b": sol["b"], "v_lsr": v_lsr,
                     "lsr_corr": sol["lsr_corr"]}
        if isinstance(self.source, SimSource):
            self.source.f_peak_baseband = f_obs - self.source.center
        self.info_label.setText(
            f"l={sol['l']:.1f}°  b={sol['b']:.1f}°\n"
            f"LST={sol['lst']}\n"
            f"Sim inject: v_lsr={v_lsr:+.1f} km/s\n"
            f"LSR corr={sol['lsr_corr']:+.2f} km/s")

        # Traffic light: |b| small = the plane is crossing your beam NOW.
        l_deg, b_deg = sol["l"], sol["b"]
        if abs(b_deg) < 2.0:
            self._set_banner("#1d7a33",
                             f"l = {l_deg:.1f}°   b = {b_deg:+.2f}°\n"
                             "ON PLANE — log / record now")
        elif abs(b_deg) < 5.0:
            self._set_banner("#9a6b00",
                             f"l = {l_deg:.1f}°   b = {b_deg:+.2f}°\n"
                             "Near plane — transit approaching")
        else:
            self._set_banner("#8a2727",
                             f"l = {l_deg:.1f}°   b = {b_deg:+.2f}°\n"
                             "Off plane — wait for transit")

    # ------------------------------------------------ target planner -------
    def plan_target(self):
        """Tell the user which altitude & local time hits galactic longitude l."""
        try:
            p = plan_pointing(self.target_l_spin.value(),
                              self.lat_spin.value(), self.lon_spin.value())
        except Exception as exc:
            self.plan_label.setText(f"Planner failed: {exc}")
            return
        self._planned_alt = p["alt_transit"] if p["observable"] else None
        l = self.target_l_spin.value()
        if not p["transits_south"]:
            self.plan_label.setText(
                f"l={l:.0f}° (Dec={p['dec']:.1f}°) transits NORTH of your zenith "
                f"— not reachable with a due-south beam from lat "
                f"{self.lat_spin.value():.1f}°. Try a lower longitude.")
            return
        if not p["observable"]:
            self.plan_label.setText(
                f"l={l:.0f}° (Dec={p['dec']:.1f}°) never rises above the southern "
                f"horizon from your latitude.")
            return
        self.plan_label.setText(
            f"l={l:.0f}° → set Altitude ≈ {p['alt_transit']:.1f}° (due south).\n"
            f"Dec={p['dec']:.1f}°.  Crosses the meridian in "
            f"{p['hours_to_transit']:.1f} h, at "
            f"{p['transit_local']:%a %H:%M} local.\n"
            f"Click “Apply altitude” to load it.")

    def apply_planned_altitude(self):
        """Always recompute for the l currently in the box — never apply a
        stale altitude from an earlier "Find pointing" click."""
        l = self.target_l_spin.value()
        try:
            p = plan_pointing(l, self.lat_spin.value(), self.lon_spin.value())
        except Exception as exc:
            self.plan_label.setText(f"Planner failed: {exc}")
            return
        if not p["observable"]:
            self._planned_alt = None
            self.plan_label.setText(
                f"l={l:.0f}° is not reachable due south from your latitude — "
                f"run “Find pointing” for details.")
            return
        alt = p["alt_transit"]
        self._planned_alt = alt
        if not self.alt_spin.isEnabled():
            self._set_status("Altitude locked while integrating — "
                             "STOP & RECORD first.")
            return
        self.alt_spin.setValue(alt)
        self._set_status(f"Altitude set to {alt:.1f}° for l={l:.0f}° "
                         f"(transits {p['transit_local']:%a %H:%M}).")

    def whats_up_tonight(self):
        """One-click list of the longitudes transiting due south in the next
        12 h, with the altitude to set and the expected Doppler sign."""
        lat, lon = self.lat_spin.value(), self.lon_spin.value()
        rows = []
        for l in range(0, 360, 10):
            try:
                p = plan_pointing(float(l), lat, lon)
            except Exception:
                continue
            if not p["observable"] or p["hours_to_transit"] > 12.0:
                continue
            if l < 90:
                tag = "redshifted (inner galaxy)"
            elif 165 <= l <= 195:
                tag = "v ≈ 0 (anticenter)"
            else:
                tag = "outer galaxy (mild)"
            rows.append((p["hours_to_transit"],
                         f"l={l:3d}°   alt {p['alt_transit']:5.1f}°   "
                         f"{p['transit_local']:%a %H:%M} (+{p['hours_to_transit']:4.1f} h)   "
                         f"{tag}"))
        rows.sort()
        text = ("Set the dish to 'alt', wait for the listed time, then log "
                "when the banner turns green:\n\n"
                + "\n".join(r[1] for r in rows)) if rows else \
            "Nothing on the plane transits due south in the next 12 h."

        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle("What can I see tonight?")
        dlg.resize(560, 420)
        v = QtWidgets.QVBoxLayout(dlg)
        box = QtWidgets.QPlainTextEdit(text)
        box.setReadOnly(True)
        box.setStyleSheet("font-family: Consolas, monospace; font-size: 12px;")
        v.addWidget(box)
        btn = QtWidgets.QPushButton("Close")
        btn.clicked.connect(dlg.accept)
        v.addWidget(btn)
        dlg.exec()

    # -------------------------------------------------------- actions ------
    def start_preview(self):
        # stop any existing stream first
        self._stop_worker()
        kind = self.source_combo.currentText()
        gain = self.gain_combo.currentData()
        try:
            self.source = make_source(kind, gain)
        except Exception as exc:
            self._set_status(f"Could not create {kind} source: {exc}")
            return
        self._on_sim_changed()
        self.update_injection()

        self.worker = AcquisitionWorker(self.source)
        self.worker.spectrum_ready.connect(self.on_spectrum)
        self.worker.status.connect(self._set_status)
        self.worker.failed.connect(self._on_worker_failed)
        self.worker.start()

        self.btn_log.setEnabled(True)
        self.btn_record.setEnabled(False)
        self.btn_baseline.setEnabled(True)
        self._set_status(f"Live preview started ({kind}).")

    def start_logging(self):
        if self.worker is None:
            self._set_status("Start the live preview first.")
            return
        self.is_logging = True
        self.log_sum = None
        self.log_count = 0
        self.integration_count = 0
        self._log_t0 = time.time()
        self._log_last_data_t = None
        self.btn_record.setEnabled(True)
        self.btn_log.setEnabled(False)
        # Changing the pointing mid-integration would silently corrupt the scan.
        self.alt_spin.setEnabled(False)
        self.integ_label.setText("Integrating: 0 s  (0 blocks)")
        self._set_status("LOGGING drift scan — integrating spectra…")

    def stop_and_record(self):
        # Integration time = span of ACTUALLY RECEIVED data, not wall clock at
        # the moment the button is pressed.  If the source died mid-scan (SDR
        # unplugged, USB fault) the clock would otherwise keep "integrating"
        # nothing — the July 15 2026 19.4-hour phantom scan.
        t0 = self._log_t0
        t_last = self._log_last_data_t if self._log_last_data_t else t0
        integ_s = max(0.0, t_last - t0) if t0 else 0.0
        # Coordinates/timestamp belong to the MIDDLE of the data span (the sky
        # drifts ~0.25 deg/min), not to whenever the row happens to be written.
        t_mid = (t0 + 0.5 * integ_s) if t0 else None
        self.is_logging = False
        self._log_t0 = None
        self._log_last_data_t = None
        self.btn_log.setEnabled(True)
        self.btn_record.setEnabled(False)
        self.alt_spin.setEnabled(True)

        if self.log_count == 0 or self.log_sum is None:
            self._set_status("No integrated data to record.")
            return

        # Final integration: average LINEAR power (radiometry), dB for humans.
        avg_lin = self.log_sum / self.log_count
        avg_db = 10.0 * np.log10(np.maximum(avg_lin, 1e-30))
        if self.chk_subtract.isChecked() and self.baseline_lin is not None \
                and len(self.baseline_lin) == len(avg_lin):
            # (ON - OFF) / OFF: flat, dimensionless line-to-continuum ratio.
            work = (avg_lin - self.baseline_lin) / np.maximum(self.baseline_lin,
                                                              1e-30)
        else:
            work = avg_lin

        # ---- restrict to the search window ----
        freqs = self.freqs
        mask = np.abs(freqs - F_REST_HZ) <= SEARCH_HALF_HZ
        if not np.any(mask):
            self._set_status("Search window empty — check frequency axis.")
            return
        idx_in = np.where(mask)[0]
        bin_hz = float(freqs[1] - freqs[0])

        # ---- smooth to HI line width before detecting ----
        # A real HI line is 50-500 kHz wide, i.e. spread over ~100-900 of these
        # 586 Hz bins; per-bin it is invisible under radiometer noise.  Boxcar-
        # smooth to SMOOTH_DETECT_HZ so broad-but-faint beats narrow-but-noisy.
        k = max(1, int(round(SMOOTH_DETECT_HZ / bin_hz)))
        if k % 2 == 0:
            k += 1
        kern = np.ones(k) / k
        work_sm = np.convolve(np.pad(work, k // 2, mode="edge"), kern,
                              mode="valid")

        win_find = work_sm[idx_in]                   # smoothed, for detection
        win_raw = work[idx_in]                       # per-bin, for the RFI veto
        win_db = avg_db[idx_in]                      # true power, for reporting

        # ---- multi-peak detection: one component per arm on the sightline ----
        med = float(np.median(win_find))
        sigma = float(np.median(np.abs(win_find - med))) * 1.4826 + 1e-30
        min_sep = max(5, int(round(40e3 / bin_hz)))  # >=~40 kHz apart (arms)
        peaks, _props = find_peaks(win_find, prominence=5.0 * sigma,
                                   distance=min_sep, width=3)

        # Narrowband-RFI veto: a real line spreads its power across the whole
        # smoothing kernel; a coherent spike concentrates it in ~1 raw bin.
        def _is_narrow_rfi(p: int) -> bool:
            lo, hi = max(0, p - k // 2), min(len(win_raw), p + k // 2 + 1)
            seg = np.clip(win_raw[lo:hi] - med, 0.0, None)
            tot = float(seg.sum())
            return tot > 0.0 and float(seg.max()) / tot > 0.5
        peaks = np.array([p for p in peaks if not _is_narrow_rfi(p)],
                         dtype=int)

        if len(peaks) == 0:                          # always record the best bin
            peaks = np.array([int(np.argmax(win_find))])
        order = np.argsort(win_db[peaks])[::-1]      # strongest first
        peaks = peaks[order][:MAX_COMPONENTS]

        try:                                         # FWHM of each retained peak
            w_meas_hz = peak_widths(win_find, peaks, rel_height=0.5)[0] * bin_hz
            # remove the boxcar's own broadening (quadrature approximation)
            fwhm_khz = np.sqrt(np.maximum(w_meas_hz ** 2 - (k * bin_hz) ** 2,
                                          0.0)) / 1e3
        except Exception:
            fwhm_khz = np.zeros(len(peaks))

        # ---- coordinates (shared by all components of this scan) ----
        # Solved at the data-span midpoint (self-test overrides via _fixed_time).
        when = self._fixed_time
        if when is None and t_mid is not None:
            when = Time(t_mid, format="unix", scale="utc")
        try:
            sol = astro_solve(self.lat_spin.value(), self.lon_spin.value(),
                              self.alt_spin.value(), when=when)
        except Exception as exc:
            self._set_status(f"Astro solve failed at record: {exc}")
            return

        src_name = "Simulation" if isinstance(self.source, SimSource) else "Hardware"

        # ---- save the full averaged spectrum (the real science product) ----
        spec_rel = self._save_spectrum(freqs, avg_lin, sol, src_name, integ_s)

        n = len(peaks)
        v_list = []
        best_snr = 0.0
        for i, p in enumerate(peaks):
            gidx = int(idx_in[p])
            # Sub-bin peak frequency: 3-point parabolic interpolation on the
            # smoothed spectrum (the raw per-bin spectrum is noise-dominated).
            delta = 0.0
            if 0 < gidx < len(work_sm) - 1:
                y0, y1, y2 = work_sm[gidx - 1], work_sm[gidx], work_sm[gidx + 1]
                denom = y0 - 2.0 * y1 + y2
                if abs(denom) > 1e-30:
                    delta = float(np.clip(0.5 * (y0 - y2) / denom, -1.0, 1.0))
            f_obs = float(freqs[gidx]) + delta * bin_hz
            power = float(avg_db[gidx])              # true power, not baselined
            snr = (float(win_find[p]) - med) / sigma
            best_snr = max(best_snr, snr)
            v_topo = C_KMS * (F_REST_HZ - f_obs) / F_REST_HZ    # radio convention
            v_lsr = v_topo + sol["lsr_corr"]
            v_list.append(v_lsr)

            row = {
                "Timestamp_UTC": sol["utc"],
                "Source": src_name,
                "Lat_GPS": f"{self.lat_spin.value():.5f}",
                "Lon_GPS": f"{self.lon_spin.value():.5f}",
                "Dish_Azimuth": "180",
                "Dish_Altitude": f"{self.alt_spin.value():.2f}",
                "Galactic_Longitude_L": f"{sol['l']:.4f}",
                "Galactic_Latitude_B": f"{sol['b']:.4f}",
                "Peak_Frequency_MHz": f"{f_obs / 1e6:.6f}",
                "Peak_Power_dB": f"{power:.3f}",
                "Doppler_Velocity_Topo_km_s": f"{v_topo:.3f}",
                "Doppler_Velocity_LSR_km_s": f"{v_lsr:.3f}",
                "Component_Index": str(i + 1),
                "N_Components": str(n),
                "SNR_sigma": f"{snr:.1f}",
                "Line_FWHM_kHz": f"{float(fwhm_khz[i]):.1f}",
                "Integration_s": f"{integ_s:.1f}",
                "Spectrum_File": spec_rel,
            }
            self._append_csv(row)
            self.records.append({"l": sol["l"], "b": sol["b"],
                                 "power": power, "v_lsr": v_lsr})
            self._add_point_to_maps(sol["l"], sol["b"], power, v_lsr)

        self._refresh_maps()

        # clear temp integration arrays
        self.log_sum = None
        self.log_count = 0
        self.integ_label.setText(
            f"Last scan: {integ_s:.0f} s, {n} component(s), "
            f"best SNR {best_snr:.0f}σ")

        vstr = ", ".join(f"{v:+.0f}" for v in v_list)
        self._set_status(
            f"Recorded {n} component(s) [{src_name}] at l={sol['l']:.1f}°, "
            f"b={sol['b']:.1f}°: v_lsr = [{vstr}] km/s")

    def _save_spectrum(self, freqs, avg_lin, sol, src_name, integ_s) -> str:
        """Write the scan's averaged linear spectrum to spectra/*.npz and return
        the app-dir-relative path ('unsaved' if the write fails)."""
        try:
            os.makedirs(SPECTRA_DIR, exist_ok=True)
            ts = "".join(ch for ch in sol["utc"] if ch.isdigit() or ch == "T")[:15]
            path = os.path.join(SPECTRA_DIR, f"scan_{ts}.npz")
            k = 2
            while os.path.exists(path):              # same-second collision
                path = os.path.join(SPECTRA_DIR, f"scan_{ts}_{k}.npz")
                k += 1
            baseline = self.baseline_lin if (
                self.baseline_lin is not None
                and len(self.baseline_lin) == len(avg_lin)) else np.array([])
            np.savez_compressed(
                path, freqs_hz=freqs, psd_lin=avg_lin, baseline_lin=baseline,
                l=sol["l"], b=sol["b"], lsr_corr=sol["lsr_corr"],
                utc=sol["utc"], source=src_name, integration_s=integ_s)
            return os.path.relpath(path, _app_dir())
        except Exception as exc:
            self._set_status(f"Spectrum save failed (row still logged): {exc}")
            return "unsaved"

    def _on_subtract_toggled(self, *_):
        # Raw view lives near -70 dB, subtracted view near 0 dB: a frozen
        # y-axis would show an apparently blank plot after toggling, so
        # force one rescale whenever the view mode changes.
        self._force_rescale = True
        self._redraw_spectrum()

    def capture_baseline(self):
        if self.ema_lin is None:
            self._set_status("No spectrum yet — start the preview first.")
            return
        # The OFF spectrum divides every later scan, so integrate a real one
        # (BASELINE_BLOCKS) instead of snapshotting one noisy display frame.
        self._bl_sum = None
        self._bl_count = 0
        self._bl_target = BASELINE_BLOCKS
        self.btn_baseline.setEnabled(False)
        self._set_status(f"Capturing cold-sky baseline — integrating "
                         f"{BASELINE_BLOCKS} blocks, hold the dish still…")

    # -------------------------------------------------- spectrum handling --
    def on_spectrum(self, freqs, power_db):
        """GUI-thread slot: all accumulation & drawing happen here (no locks)."""
        self.freqs = freqs
        self.last_spectrum = power_db
        self.spectrum_count += 1

        # All accumulation happens in LINEAR power; dB is display-only.
        p_lin = np.power(10.0, power_db / 10.0)

        # light EMA smoothing for a stable live display
        if self.ema_lin is None or len(self.ema_lin) != len(p_lin):
            self.ema_lin = p_lin.copy()
        else:
            self.ema_lin = 0.6 * self.ema_lin + 0.4 * p_lin

        # baseline capture in progress: accumulate a proper OFF integration
        if self._bl_target > 0:
            if self._bl_sum is None or len(self._bl_sum) != len(p_lin):
                self._bl_sum = p_lin.copy()
                self._bl_count = 1
            else:
                self._bl_sum += p_lin
                self._bl_count += 1
            if self._bl_count % 25 == 0:
                self._set_status(f"Capturing baseline… {self._bl_count}/"
                                 f"{self._bl_target} blocks")
            if self._bl_count >= self._bl_target:
                self.baseline_lin = self._bl_sum / self._bl_count
                self._bl_sum = None
                self._bl_target = 0
                self.btn_baseline.setEnabled(True)
                self._force_rescale = True
                self.chk_subtract.setChecked(True)
                self._set_status(
                    f"Cold-sky baseline captured ({self._bl_count} blocks). "
                    "Subtraction enabled — leave it on; the flat ~0 dB line "
                    "is correct.")

        if self.is_logging:
            if self.log_sum is None or len(self.log_sum) != len(p_lin):
                self.log_sum = p_lin.copy()
                self.log_count = 1
            else:
                self.log_sum += p_lin
                self.log_count += 1
            self.integration_count = self.log_count
            self._log_last_data_t = time.time()
            elapsed = self._log_last_data_t - self._log_t0 if self._log_t0 else 0.0
            self.integ_label.setText(
                f"Integrating: {elapsed:.0f} s  ({self.log_count} blocks)")

        self._redraw_spectrum()

    def _redraw_spectrum(self):
        if self.freqs is None or self.ema_lin is None:
            return
        # show integrated average while logging, else the smoothed live view
        if self.is_logging and self.log_count > 0:
            shown_lin = self.log_sum / self.log_count
        else:
            shown_lin = self.ema_lin
        if self.chk_subtract.isChecked() and self.baseline_lin is not None \
                and len(self.baseline_lin) == len(shown_lin):
            # ratio in dB == the familiar "baseline-subtracted" flat view
            shown = 10.0 * np.log10(np.maximum(shown_lin, 1e-30)
                                    / np.maximum(self.baseline_lin, 1e-30))
        else:
            shown = 10.0 * np.log10(np.maximum(shown_lin, 1e-30))

        x = self.freqs / 1e6
        self.spec_line.set_data(x, shown)

        # overlay smoothed to HI line width (same kernel as detection)
        if len(self.freqs) > 1:
            bin_hz = float(self.freqs[1] - self.freqs[0])
            k = max(1, int(round(SMOOTH_DETECT_HZ / bin_hz)))
            if k % 2 == 0:
                k += 1
        else:
            k = 1
        if k > 1:
            sm = np.convolve(np.pad(shown, k // 2, mode="edge"),
                             np.ones(k) / k, mode="valid")
        else:
            sm = shown
        self.spec_line_sm.set_data(x, sm)

        # mark the current peak within the search window (on the smoothed
        # curve — that is where a real line lives)
        mask = np.abs(self.freqs - F_REST_HZ) <= SEARCH_HALF_HZ
        if np.any(mask):
            idx_in = np.where(mask)[0]
            pk = idx_in[int(np.argmax(sm[idx_in]))]
            self.peak_marker.set_data([x[pk]], [sm[pk]])

        if not self.chk_freeze_y.isChecked() or self._force_rescale:
            ax = self.c_spec.ax
            ymin, ymax = float(np.min(shown)), float(np.max(shown))
            pad = max(3.0, 0.1 * (ymax - ymin))
            ax.set_ylim(ymin - pad, ymax + pad)
            self._force_rescale = False
        self.c_spec.draw_idle()

    # ----------------------------------------------------------- maps ------
    def _add_point_to_maps(self, l, b, power, v_lsr):
        li = int(l) % L_BINS
        bi = int(round(b)) - B_MIN
        if 0 <= bi < B_BINS:
            # keep the strongest detection in each cell
            cur = self.lb_map[bi, li]
            if np.isnan(cur) or power > cur:
                self.lb_map[bi, li] = power
        self.lv_l.append(l)
        self.lv_v.append(v_lsr)
        self.lv_p.append(power)

    def _refresh_maps(self):
        # 2D l-b heatmap
        self.map_im.set_data(self.lb_map)
        if np.isfinite(self.lb_map).any():
            vmin = float(np.nanmin(self.lb_map))
            vmax = float(np.nanmax(self.lb_map))
            if vmax - vmin < 1.0:
                vmax = vmin + 1.0
            self.map_im.set_clim(vmin, vmax)
            self.map_cbar.update_normal(self.map_im)
        self.c_map.draw_idle()

        # longitude-velocity scatter
        if self.lv_l:
            self.lv_scatter.set_offsets(np.column_stack([self.lv_l, self.lv_v]))
            self.lv_scatter.set_array(np.asarray(self.lv_p))
            pmin, pmax = min(self.lv_p), max(self.lv_p)
            if pmax - pmin < 1.0:
                pmax = pmin + 1.0
            self.lv_scatter.set_clim(pmin, pmax)
            self.lv_cbar.update_normal(self.lv_scatter)
        self.c_lv.draw_idle()

    # ---------------------------------------------------------- misc -------
    def _append_csv(self, row: dict):
        new_file = not os.path.exists(CSV_PATH)
        with open(CSV_PATH, "a", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=CSV_HEADER)
            if new_file:
                writer.writeheader()
            writer.writerow(row)

    def _set_status(self, msg: str):
        if hasattr(self, "status_lbl"):
            self.status_lbl.setText(msg)
        # Console echo must never crash the GUI: in a --windowed / pythonw build
        # sys.stdout is None, and a plain console may be cp1252 (no Unicode).
        try:
            if sys.stdout is not None:
                line = f"[status] {msg}"
                enc = getattr(sys.stdout, "encoding", None) or "ascii"
                sys.stdout.write(line.encode(enc, "replace").decode(enc) + "\n")
        except Exception:
            pass

    def _on_worker_failed(self, msg: str):
        self._stop_worker()
        if self._bl_target > 0:            # baseline capture died with the source
            self._bl_target = 0
            self._bl_sum = None
            self.btn_baseline.setEnabled(True)
        if self.is_logging:
            # The source died mid-scan (e.g. SDR unplugged).  Do NOT leave the
            # integration "running" against a dead source: finalize now, while
            # the pointing/time still match the data, using only the data span
            # actually received.  With nothing accumulated, just reset the UI.
            if self.log_count > 0 and self.log_sum is not None:
                self.stop_and_record()
                self._set_status(
                    f"ERROR: {msg} — source lost mid-scan; the integration was "
                    f"auto-recorded from the data received before the failure.")
            else:
                self.is_logging = False
                self._log_t0 = None
                self._log_last_data_t = None
                self.btn_log.setEnabled(True)
                self.btn_record.setEnabled(False)
                self.alt_spin.setEnabled(True)
                self.integ_label.setText("Scan aborted — source lost.")
                self._set_status(f"ERROR: {msg} — scan aborted (no data).")
        else:
            self._set_status(f"ERROR: {msg}")

    def _stop_worker(self):
        if self.worker is not None:
            try:
                self.worker.spectrum_ready.disconnect(self.on_spectrum)
            except Exception:
                pass
            self.worker.stop()
            self.worker = None
        self.source = None

    def closeEvent(self, event):
        self._stop_worker()
        super().closeEvent(event)


# ===========================================================================
#  HEADLESS SELF-TEST
# ===========================================================================
def run_self_test(app: QtWidgets.QApplication, win: MainWindow) -> int:
    """Drive the whole pipeline in Simulation mode and verify the outputs."""
    state = {"step": 0, "t0": time.time(), "result": None}

    def fail(msg):
        print(f"[self-test] FAIL: {msg}")
        state["result"] = 1
        app.quit()

    def done_ok():
        print("[self-test] PASS: live spectrum, CSV row, and both maps verified.")
        state["result"] = 0
        app.quit()

    def tick():
        # global timeout guard
        if time.time() - state["t0"] > 25:
            return fail("timed out")

        step = state["step"]
        try:
            if step == 0:
                print("[self-test] step 0: start live preview (Simulation)")
                win.source_combo.setCurrentText("Simulation")
                # Deterministic sky: fix the time and point at the galactic
                # plane (|b|<20) so the row lands on BOTH maps.
                win._fixed_time = Time("2026-07-15T04:00:00", scale="utc")
                lat, lon = win.lat_spin.value(), win.lon_spin.value()
                best_alt, best_b = 40.0, 999.0
                for alt in range(3, 88):
                    b = abs(astro_solve(lat, lon, float(alt),
                                        when=win._fixed_time)["b"])
                    if b < best_b:
                        best_b, best_alt = b, float(alt)
                win.alt_spin.setValue(best_alt)
                print(f"[self-test]   fixed UTC=2026-07-15T04:00, "
                      f"altitude={best_alt:.0f}° -> |b|={best_b:.2f}°")
                win.start_preview()
                state["step"] = 1

            elif step == 1:
                if win.spectrum_count >= 4:
                    print(f"[self-test] step 1: {win.spectrum_count} spectra received "
                          "-> live simulated spectrum is displaying")
                    win.start_logging()
                    state["step"] = 2

            elif step == 2:
                if win.integration_count >= 6:
                    print(f"[self-test] step 2: integrated {win.integration_count} "
                          "blocks -> STOP & RECORD ROW")
                    n_before = len(win.records)
                    win.stop_and_record()
                    if len(win.records) < n_before + 1:
                        return fail("no record was appended")
                    state["step"] = 3

            elif step == 3:
                # verify CSV
                if not os.path.exists(CSV_PATH):
                    return fail("CSV file was not created")
                with open(CSV_PATH, newline="") as fh:
                    rows = list(csv.DictReader(fh))
                if len(rows) < 1:
                    return fail("CSV has no data rows")
                last = rows[-1]
                if list(last.keys()) != CSV_HEADER:
                    return fail(f"CSV header mismatch: {list(last.keys())}")
                for k in CSV_HEADER:
                    if last[k] in ("", None):
                        return fail(f"CSV field '{k}' is empty")
                fmhz = float(last["Peak_Frequency_MHz"])
                if not (1419.0 <= fmhz <= 1421.5):
                    return fail(f"peak freq {fmhz} MHz outside search window")
                if float(last["SNR_sigma"]) < 3.0:
                    return fail(f"detection SNR too low: {last['SNR_sigma']}σ")
                if float(last["Integration_s"]) <= 0.0:
                    return fail("Integration_s not positive")

                # verify the saved spectrum file
                spath = os.path.join(_app_dir(), last["Spectrum_File"])
                if not os.path.exists(spath):
                    return fail(f"spectrum file missing: {spath}")
                with np.load(spath) as d:
                    if "psd_lin" not in d.files or d["psd_lin"].size == 0:
                        return fail("spectrum npz has no psd_lin data")

                # verify both maps got a point
                if not np.isfinite(win.lb_map).any():
                    return fail("l-b heatmap has no finite cell")
                if len(win.lv_l) < 1:
                    return fail("longitude-velocity diagram has no point")

                print("[self-test] recorded row:")
                for k in CSV_HEADER:
                    print(f"           {k:28s} = {last[k]}")
                print(f"[self-test] l-b map finite cells : {int(np.isfinite(win.lb_map).sum())}")
                print(f"[self-test] l-v scatter points   : {len(win.lv_l)}")
                done_ok()

        except Exception:
            traceback.print_exc()
            return fail("exception during self-test")

    timer = QTimer()
    timer.timeout.connect(tick)
    timer.start(120)

    app.exec()
    return state["result"] if state["result"] is not None else 1


# ===========================================================================
#  ENTRY POINT
# ===========================================================================
def main():
    global CSV_PATH, SPECTRA_DIR
    self_test = "--self-test" in sys.argv
    force_sim = self_test or ("--sim" in sys.argv)

    if self_test:
        # No display needed for the automated verification.
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        # Never pollute the real observing log with self-test rows: redirect
        # both outputs to disposable selftest_* paths and start them clean.
        CSV_PATH = os.path.join(_app_dir(), "selftest_log.csv")
        SPECTRA_DIR = os.path.join(_app_dir(), "selftest_spectra")
        for p in (CSV_PATH, CSV_PATH + ".bak"):
            if os.path.exists(p):
                os.remove(p)
        shutil.rmtree(SPECTRA_DIR, ignore_errors=True)

    app = QtWidgets.QApplication(sys.argv)
    apply_dark_theme(app)
    win = MainWindow(force_source="sim" if force_sim else None)
    win.show()

    if self_test:
        code = run_self_test(app, win)
        win._stop_worker()
        sys.exit(code)
    else:
        sys.exit(app.exec())


if __name__ == "__main__":
    main()
