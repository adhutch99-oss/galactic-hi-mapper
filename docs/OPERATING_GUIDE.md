# Galactic Plane 21-cm HI Mapper — Complete Operating & Context Guide

> **Purpose of this file:** the complete operating manual — hardware, observing
> procedure, science background, data schema and code architecture — written so
> someone with no prior context can build the instrument, run a session and
> understand every number the software produces.
>
> Observing-site coordinates in this document are rounded to 0.1° (~10 km) for
> privacy. Only the latitude affects the pointing geometry; the longitude only
> shifts transit clock times.

---

## 1. What this project is

A single-window desktop app for a DIY radio telescope that detects and maps the
**21-cm neutral hydrogen line** (rest frequency **1420.405752 MHz**) from the
galactic plane, to reveal **galactic rotation and spiral-arm structure**.

Built and operated by one person as an amateur astronomy project. The software was
developed and validated in **Simulation mode** before the hardware arrived, so
switching to the real receiver is a single dropdown change; the simulator is still
used for regression testing (`--self-test`).

---

## 2. Files (all in the repository root)

| File | What it is |
|---|---|
| `hydrogen_mapper.py` | The whole app — single self-contained file (PyQt6 GUI, dark theme). |
| `analyze_log.py` | Companion: CSV (+ saved spectra) → l–v scatter, 2-D l–v image, rotation curve, face-on galaxy map. |
| `galactic_plane_log.csv` | The observing log. One row per detected component. Auto-created/loaded. |
| `spectra/*.npz` | Full averaged spectrum of every recorded scan (auto-saved; feeds the 2-D l–v image). |
| `Launch Hydrogen Mapper.bat` | Double-click to open the app (runs the correct Python). |
| `Make Plots (Analyze Log).bat` | Double-click to generate + show the analysis plots (all rows). |
| `Make Plots (Hardware Only).bat` | Same, but ignores `Source=Simulation` rows — use once real hardware is running. |
| `OPERATING_GUIDE.md` | This file. |
| `GRAPHS_GUIDE.txt` | Plain-text walkthrough: how to make the plots, how to read each one, green-line vs green-banner explainer, resume checklist. |
| `selftest_log.csv`, `selftest_spectra/` | Disposable outputs of `--self-test` (it never touches the real log). |

---

## 3. How to run it

**Interpreter (important):** use **Python 3.13**. At the time of writing PyQt6 /
scipy / astropy wheels were not yet published for 3.14, so a 3.14 interpreter
will fail on import. `pip install -r requirements.txt` into a 3.13 environment.

- **Normal use:** double-click **`Launch Hydrogen Mapper.bat`**.
- **Manual:** `python hydrogen_mapper.py`
- **Headless self-test (verifies the whole pipeline):**
  `python313 hydrogen_mapper.py --self-test` → prints `PASS`, exits 0. It runs
  offscreen in Simulation, records a row, and checks the CSV row (incl. SNR /
  FWHM / Integration_s), the saved spectrum `.npz`, and both maps. It pins a
  fixed UTC (2026-07-15T04:00) and sweeps altitude to hit the plane so it's
  deterministic. It writes to `selftest_log.csv` / `selftest_spectra/` — the
  real observing log is never touched by testing.

**Installed into Python 3.13:** numpy, scipy, matplotlib, PyQt6, astropy,
**pyrtlsdr**, **pyrtlsdrlib** (bundles the `librtlsdr` DLL). So when hardware arrives,
**no more Python installs are needed.**

---

## 4. Hardware (the real backend)

Signal chain, in order:
`Dish feed → LMR400 cable → SAWbird+ H1 LNA → RG316 pigtail → NESDR SMArTee v2 SDR → USB`.

- SDR: **Nooelec NESDR SMArTee v2** (RTL2832U + R820T2). Its SMA port has an
  **ALWAYS-ON 4.5 V bias-tee** — there is **nothing to switch on**, and the code
  deliberately never toggles a bias-tee.
- LNA: **Nooelec SAWbird+ H1** (+40 dB, SAW filter @ 1420 MHz), powered by that
  bias-tee. **Must stay inline.**
- Antenna: **Nooelec 20 dBi parabolic mesh dish** (1.4 GHz). Beam ~8–10° wide.

**The ONE remaining setup step when parts arrive** (can't be done without the device):
install the Windows **WinUSB** driver via **Zadig** — plug in the SDR, run Zadig,
Options → List All Devices → select the RTL "Bulk-In, Interface (Interface 0)" →
set target driver to **WinUSB** → Replace/Install. One-time per computer.

If no SDR is found at startup the app **auto-falls back to Simulation** and says so.

---

## 5. The instrument concept (READ THIS — it's the source of most confusion)

This is a **meridian drift-scan** telescope:

- The dish points **due south (Azimuth = 180°, fixed)**. You never steer in azimuth.
- The **only aiming control is Altitude** (elevation angle, read off a protractor).
- **Altitude picks a horizontal stripe of sky (a declination). The current TIME picks
  which galactic longitude of that stripe is sitting due south.** You need BOTH right
  to hit a target — altitude alone is not enough.
- As Earth rotates, different galactic longitudes drift across the southern meridian
  (~10–15° of galactic longitude per hour). You catch each one as it transits.

**`b` (galactic latitude) = how far your beam is above/below the galactic plane.**
- `b ≈ 0` means you're centered ON the plane (where the gas/arms are) → strong signal.
- `b` far from 0 = off the plane → weak/no signal.
- For any target longitude, `b` sweeps through 0 exactly at that longitude's **transit
  time**. So "wait for b≈0" and "record at transit time" are the same instant.
- **b has nothing to do with redshift/blueshift.** That's a common user misconception.
  b = pointing; redshift = velocity. Separate things.

**Fixed site (Long Island):** observations are all taken from **lat 40.7° N (approx.)
lon −73.6° W (approx.)** (West = negative). These are baked in as the app defaults.
Because latitude is fixed, **each galactic longitude maps to ONE fixed altitude
forever** — only the transit *time* shifts ~4 min earlier each night.

### Longitude → dish-altitude table (for lat 40.7° N, due south)

Reachable when target declination ≤ 40.7°. Inner galaxy (rotation-curve gold):

| l (deg) | Dish altitude | l (deg) | Dish altitude |
|---:|---:|---:|---:|
| 0 | 20° (low, near horizon) | 50 | 64° |
| 10 | 29° | 60 | 73° |
| 20 | 38° | 70 | 82° |
| 30 | 47° | 80 | 90° (zenith) |
| 40 | 56° | | |

Outer galaxy / anticenter (reachable, low-velocity): l=170→87°, 180→78°, 200→61°,
220→43°, 240→25°, 260→9°, 270→1° (horizon).

**Unreachable from Long Island:** l ≈ 90°–165° (transits NORTH of zenith — Cygnus/
Cassiopeia/Perseus, dec too high) and l ≈ 275°–335° (below the southern horizon).

So the realistic mapping range is **inner galaxy l = 10°–80°** (all redshifted) plus
the anticenter/outer l = 170°–265°. Inner galaxy transits ~1–2 AM in summer — high and easy.

### Setting the dish altitude (physical measurement)

The app's "Altitude (deg)" is the **beam elevation**. Measure it along the **feed
boom** (the arm holding the LNA) — that arm points down the beam (0°=horizon, 90°=up).
Do NOT measure a random spot on the rim.

1. **Easiest — phone inclinometer/level app** laid flat on the feed boom (or a mount
   surface parallel to it). Reads degrees directly, ~1° accuracy. Good enough (beam is ~8–10° wide).
2. **Tape-measure / trig** — run `angle_calc.py` (or `Measure Altitude Angle.bat`):
   - *Board & drop (easiest):* rest a straight level of known **length** along the boom,
     measure the vertical **height difference** between its ends → angle = asin(drop/length).
   - *Rise & run:* angle = atan(rise/run).
3. **Calibrate with the Sun (definitive):** the Sun is a strong 1420 MHz source. With it
   due south, sweep altitude to maximize signal; that altitude should equal the Sun's
   known elevation at that timestamp (from a planetarium app). Any consistent difference
   is your boom-vs-beam offset — note it and apply it to every pointing.

Tools: `angle_calc.py`, `Measure Altitude Angle.bat`.

---

## 6. The science (what the graphs mean)

- **21-cm line:** neutral hydrogen emits at 1420.405752 MHz. Its observed frequency
  is Doppler-shifted by the gas's line-of-sight velocity.
- **Doppler / velocity:** `v = c·(f_rest − f_obs)/f_rest`, radio convention
  (**positive = receding = redshift**, negative = approaching = blueshift).
- **v_LSR:** the topocentric velocity corrected to the Local Standard of Rest
  (removes Earth's motion + the Sun's peculiar motion). This is the scientifically
  meaningful number — the app logs both topo and LSR.
- **Red vs blue = which galactic quadrant you look at, NOT which arm:**
  - Inner galaxy **Quadrant I (l ≈ 0–90°)** → **redshifted** (positive v_LSR).
  - Anticenter (l ≈ 180°) → ~0.
  - Its blueshifted mirror is Quadrant IV (l ≈ 270–360°), which sits too far south to
    see from New York.
- **Distinguishing arms:** along one sightline, each arm shows up as a **separate
  velocity peak** (different distance → different rotation velocity). Across many
  longitudes, each arm traces a continuous **lane** in the longitude–velocity diagram.
  (Caveat: kinematic distance ambiguity — one velocity = two distances except at the
  tangent point.)
- **Tangent-point method (rotation curve):** works only for the **inner galaxy**
  (0<l<90 or 270<l<360). The most extreme velocity at each longitude comes from the
  tangent point, where `R = R0·|sin l|` and `V(R) = |v_tangent| + V0·|sin l|`.
  Constants: **R0 = 8.2 kpc, V0 = 220 km/s.**

### The three live graphs (right side of the app)
1. **Top — Live spectrum:** Frequency (MHz) vs Power (dB), with a **velocity
   axis (km/s) along the top** so Doppler can be read directly. Faint blue
   line = raw per-bin spectrum; **bold gold line = the same spectrum smoothed
   to HI line width (50 kHz) — watch THIS one: a real hydrogen line appears
   here as a broad bump long before anything is visible per-bin.** **Green
   dashed line = rest frequency (1420.4058 MHz) = zero-velocity marker.** Peak
   **left** of green = redshift (receding); **right** = blueshift. Orange band
   = ±1 MHz search window; red triangle = peak of the smoothed curve. A
   "Freeze spectrum y-axis" checkbox stops the autoscale jumping.
2. **Middle — l–b heatmap:** galactic longitude (x) vs latitude −20..+20 (y),
   color = peak power. Each recorded plane point lights a cell. Gray bands mark
   longitudes unreachable from the fixed site.
3. **Bottom — Longitude–Velocity diagram:** longitude (x) vs v_LSR (y), colored by
   power. **The money plot** — arm lanes + rotation signature. Same gray
   unreachable bands.

A large color-coded **ON PLANE banner** at the top of the control panel shows
live l / b: green (|b|<2°) = the plane is in the beam, log/record now; amber
(|b|<5°) = transit approaching; red = off plane, wait.

---

## 7. Simulation vs Hardware (what's real)

- **Simulation** generates synthetic IQ: realistic noise floor, a modeled receiver
  bandpass ("gain skirt"), a residual DC spike, and **one injected HI line** whose
  velocity follows a toy flat-rotation model `v_lsr = 220·sin(l)·cos(b)`. Adjustable
  **Injected SNR** (sim-only; ignore on hardware) and optional fake RFI spikes.
- **KEY:** because the sim injects only ONE line on a smooth toy curve, the sim l–v
  diagram will only ever show a clean single-point curve — **it can NOT show real
  spiral arms.** Real structure requires real hardware. The sim is a "flight simulator"
  for the workflow and the pipeline, nothing more.
- Simulation data IS written to the same CSV and maps as real data, but every row is
  tagged in the **`Source`** column (`Simulation`/`Hardware`) so they're never confused.

---

## 8. The observing workflow (button-by-button)

**Rule of thumb:** a dot = ~3 min of button-pressing. You are not waiting on the app —
you're waiting on the **sky** to bring the next longitude onto the plane.

**Once per session (cold-sky baseline — NOT a recorded measurement):**
1. Point at empty, high, off-plane sky. **In summer, capture the baseline EARLY
   (~10–10:30 PM), pointed straight up (zenith)**: that's b ≈ +30…+40°, the
   best cold sky reachable. The later it gets, the closer zenith drifts to the
   plane (by ~1:30 AM the plane IS at zenith); for a post-midnight refresh use
   **altitude 40°** instead (b ≈ −15…−30°). Anything |b| > 15° works; > 30° is
   ideal. (This pattern shifts ~4 min earlier per night, like everything else.)
2. **START LIVE PREVIEW**, wait ~5 s for the line to settle.
3. **Capture Cold Sky Baseline** and **hold the dish still ~20–30 s**: since
   July 2026 the button integrates a real OFF measurement (`BASELINE_BLOCKS` =
   150 blocks, progress in the status bar) instead of snapshotting one noisy
   display frame — the old snapshot's noise was the dominant noise in every
   subtracted spectrum. Subtraction auto-enables when it finishes.
   *(Do NOT Log Drift Scan or Stop & Record for the baseline — that would log junk.)*
   **Recapture the baseline every ~30 min** (repoint off-plane, recapture, repoint
   back): on the July 2026 first-light session the receiver gain drifted so much
   over ~40 min that the last scan's continuum ratio sat at −0.4 instead of ~0,
   which resurrects bandpass ripple the subtraction should have cancelled.

**For each target longitude:**
4. Aim due south at that longitude's altitude (from the table / Find pointing /
   the **"What can I see tonight?"** button, which lists everything transiting
   in the next 12 h with altitudes and times).
5. Watch the banner; when it turns **green (`b ≈ 0`**, its transit), click
   **LOG DRIFT SCAN**. The altitude field locks while integrating (changing
   pointing mid-scan would corrupt it) and a counter shows elapsed seconds.
6. Let it integrate **1–3 minutes** (longer = cleaner for faint targets).
7. **STOP && RECORD ROW** → logs a row per detected component (with its SNR and
   line width), saves the full spectrum to `spectra/`, updates both maps.
8. **Nudge the Altitude up** toward the next longitude, wait for the **green
   banner** again, then Log → Record. Repeat. (No need to press START LIVE
   PREVIEW again — it keeps streaming.)

**Cadence:** hands-on ~3 min/dot; for well-separated longitudes (~5–10° apart, good
for a rotation curve) expect roughly **one dot every 20–30 min** — that's the sky
drifting, not idle time. A 2–3 hr session → ~5–8 dots spanning ~30–50° of longitude.
On real hardware, one scan can log **several** dots at once (one per arm).

**Multi-night:** the same longitude transits ~4 min earlier each night; altitudes never
change. Use **Find pointing** each night for that night's exact transit times.

### The Target Planner ("Find pointing")
It's an **advisor/calculator, not a steering wheel.** Type a galactic longitude → it
tells you (a) the dish **altitude**, (b) the local **transit time** ("in X h, at HH:MM"),
and (c) whether it's reachable due south. **Apply altitude** loads it into the field —
it **recomputes from the l currently typed in the box on every click** (an earlier
version applied a stale value from the last "Find pointing"; fixed July 2026).
It does NOT change what you're currently observing — the live `l` readout (altitude +
current time) is what actually gets logged.

---

## 9. CSV schema (`galactic_plane_log.csv`)

Columns (extended past the original spec, with automatic in-place migration + `.bak`
backup if an older-format file is found):

| Column | Meaning |
|---|---|
| `Timestamp_UTC` | Auto-stamped **midpoint of the integration's data span** (you never type this). |
| `Source` | `Simulation` or `Hardware`. |
| `Lat_GPS`, `Lon_GPS` | Observer location. |
| `Dish_Azimuth` | Always 180 (south). |
| `Dish_Altitude` | Elevation you set. |
| `Galactic_Longitude_L`, `Galactic_Latitude_B` | Where on the galaxy you looked. |
| `Peak_Frequency_MHz` | Measured line frequency. |
| `Peak_Power_dB` | Signal strength. |
| `Doppler_Velocity_Topo_km_s` | Raw velocity (reference only). |
| `Doppler_Velocity_LSR_km_s` | **The science velocity** — use this. |
| `Component_Index`, `N_Components` | Which arm/component of a multi-peak scan, and how many total. |
| `SNR_sigma` | Detection significance (peak over robust noise on the 50 kHz-smoothed spectrum, in σ). The fallback always records a row even with no detection — **filter on this** (analysis drops <5σ by default). |
| `Line_FWHM_kHz` | Measured line width at half prominence, corrected for detection smoothing. **≈0 means the fallback landed on a noise spike / bandpass artifact, not a line** — analysis drops FWHM < 20 kHz. |
| `Integration_s` | Span of **actually received data** (start of scan → last block received), not wall clock at the moment the row is written. |
| `Spectrum_File` | Relative path of the scan's saved full spectrum (`spectra/*.npz`). |

The app **reloads this file on startup** and repopulates both maps, so progress
persists across launches. To read it: open in Excel; for the l–v diagram make an XY
scatter of `Galactic_Longitude_L` vs `Doppler_Velocity_LSR_km_s`.

---

## 10. Analysis (`analyze_log.py` / `Make Plots ....bat`)

Reads the CSV (plus `spectra/*.npz`) and writes next to it:
- **`lv_diagram.png`** — longitude–velocity peak scatter with the flat-rotation
  tangent envelope overlaid (arm lanes).
- **`lv_diagram_2d.png`** — TRUE l–v intensity image rebuilt from the saved
  spectra (full line profiles, not just peaks). Skipped if no spectra exist yet.
- **`rotation_curve.png`** — tangent-point rotation curve (needs inner-galaxy points,
  0<l<90 or 270<l<360; will be **empty** if you only logged outer-galaxy longitudes).
- **`faceon_map.png`** — top-down Milky Way map from kinematic distances
  (near solution = filled, far = open, tangent points = squares; the near/far
  ambiguity is inherent to the method).

Flags: `--show` (pop up on screen), `--hardware-only` (ignore `Source=Simulation`
rows), `--min-snr N` (detection cut in σ, default 5; legacy rows without an SNR
always pass), `--fits` (also export `galactic_plane_log.fits`). Rows with
`Line_FWHM_kHz` < 20 (unphysically narrow — fallback/artifact rows) are always
dropped alongside the SNR cut.

---

## 11. Cleared-up misconceptions (every one of these bit me in practice)

- The planner is advice, not aiming. Recorded `l` = altitude + current time.
- Setting altitude alone doesn't give a target longitude — you also need the right time.
- `b ≈ 0` = on the plane (pointing), unrelated to redshift.
- The green line = zero-velocity reference; on-target inner-galaxy peaks land *left* (red) of it.
- **You do NOT wait for the spike to reach the green line.** You wait for the green
  *banner* (b≈0 = transit) to start the drift scan. The spike stays offset from the
  green line — that offset IS the Doppler measurement. A spike ON the green line
  would mean zero velocity (local foreground gas), not a successful detection.
- Baseline capture is NOT a recorded row (no Log/Record for it).
- Sim shows only ONE line on a toy curve — no real arms until hardware.
- You don't wait idle; you wait for the sky. ~3 min hands-on per dot.

---

## 12. Code architecture (for a future dev chat)

- **Constants:** `F_REST_HZ`, `CENTER_HZ`=1420.00 MHz (offset tune so the line sits
  ~0.4 MHz off the DC spike), `SAMPLE_RATE_HZ`=2.4 MSps, `CHUNK`=131072,
  `NPERSEG`=4096, `SEARCH_HALF_HZ`=1.0 MHz, `MAX_COMPONENTS`=5, gains list w/ default 40.2 dB.
- **DSP:** `compute_psd` (Welch-style averaged FFT → dB), `notch_dc` (interpolates the
  RTL DC spike). **All integration/averaging happens in LINEAR power** (dB is
  display-only): the EMA display buffer, the drift-scan accumulator, and the
  cold-sky baseline are linear PSDs — averaging dB is biased and loses SNR.
  Baseline subtraction is (ON−OFF)/OFF (line-to-continuum ratio).
  **Detection runs on a boxcar-smoothed spectrum** (`SMOOTH_DETECT_HZ` = 50 kHz,
  ~85 bins): real HI lines are 100–500 kHz wide, so per-bin they drown in
  radiometer noise — smoothing to line width is what made the July 2026
  first-light detection visible. `find_peaks` (scipy) then does multi-peak
  detection with a robust MAD-based prominence threshold, plus a
  **narrowband-RFI veto** (a candidate is rejected if one raw bin holds >50 % of
  the kernel window's power — a coherent spike, not thermal emission); falls
  back to argmax so a row is always recorded, but every row carries `SNR_sigma`
  and `Line_FWHM_kHz` so analysis can drop non-detections. Peak frequency is
  refined by 3-point parabolic interpolation on the smoothed spectrum (sub-bin);
  FWHM comes from `scipy.signal.peak_widths` with the boxcar's broadening
  removed in quadrature.
  `_save_spectrum` writes each scan's averaged linear spectrum to `spectra/*.npz`
  (keys: freqs_hz, psd_lin, baseline_lin, l, b, lsr_corr, utc, source,
  integration_s).
- **Data sources:** `DataSource` base → `SimSource` (synthetic IQ, freq-domain
  synthesis) and `RealSource` (pyrtlsdr `read_samples_async` into a bounded queue).
  `make_source()`, `hardware_available()`.
- **Threading:** `AcquisitionWorker(QThread)` runs read+DSP off the GUI thread and
  emits `spectrum_ready` (queued signal). **All accumulation, plotting, widget access
  happen on the GUI thread — no locks.** Do not touch Qt from the worker.
- **Astro:** `astro_solve(lat,lon,alt,az=180,when)` → AltAz→ICRS→galactic (l,b) plus
  LSR correction (astropy barycentric + solar motion U,V,W = 11.1, 12.24, 7.25).
  `plan_pointing(l,lat,lon)` for the planner. `toy_injected_v_lsr` for the sim.
- **MainWindow:** builds the scrollable left panel + 3 canvases; `stop_and_record`
  does the multi-peak detection and writes rows; `_ensure_csv_schema` migrates old
  logs (the pre-SNR-column real log will auto-migrate with a `.bak` on next
  launch); `update_injection` keeps the sim line tracking pointing/time and
  drives the ON-PLANE banner; `whats_up_tonight` lists the next-12 h transits;
  `_fixed_time` overrides the clock (self-test only).
- **UI/theme:** dark Fusion palette + Segoe UI 10 via `apply_dark_theme`;
  matplotlib canvases match via `style_axis`/`style_colorbar`; the spectrum has
  a secondary top axis in km/s (`freq_mhz_to_vlsr_axis`). The altitude spinner
  locks during integration.
- **Robustness:** `_set_status` is hardened against `sys.stdout is None` (pythonw /
  `--windowed` PyInstaller builds) and non-Unicode consoles.
  **Source-failure handling (July 2026, after the 19.4-hour phantom scan):** if
  the SDR dies mid-scan (unplugged, USB fault), `read_iq()` times out in ~5 s,
  the worker emits `failed`, and `_on_worker_failed` **finalizes the scan
  immediately** (auto Stop & Record with the data received so far) or aborts it
  if nothing was accumulated — the integration can no longer "run" against a
  dead source. `Integration_s` is the span from scan start to the **last block
  actually received**, and coordinates/timestamp are solved at the **midpoint of
  that data span**, so a row recorded late still describes where the dish was
  actually looking.
- **Packaging:** PyInstaller `--windowed`, collect astropy / astropy_iers_data /
  pyerfa data; hidden-import scipy.signal, astropy.coordinates, (rtlsdr if bundling HW).

---

## 13. Open ideas / not yet built

(Built in July 2026: "What can I see tonight?" button, velocity axis on the
spectrum, face-on arm map, linear-power integration, SNR/FWHM columns,
per-scan spectrum saving + 2-D l–v image, dark theme, ON-PLANE banner.)

- Near/far kinematic-distance **resolution** (the face-on map currently shows
  both solutions; resolving needs line-width/latitude-extent heuristics).
- A per-session free-text notes column.
- A PPM frequency-offset field (the SMArTee v2's 0.5 ppm TCXO makes this ~0.1 km/s,
  so it's low priority).

**Status as of July 17 2026 (see the top-level README for the current data set):** hardware arrived and works. **First-light
session July 15 2026 (~11:30 PM–12:15 AM local): four ~3-min drift scans at
l = 11°, 22°, 28°, 33°, all |b| < 1.1°.** The smoothed spectra
(`first_light_smoothed.png`) show a broad redshifted emission bump rising at
v_LSR ≈ 0 with a tail to ~+80–100 km/s in every on-plane scan, weakening with
l — consistent with real Quadrant-I galactic HI, though at ~3–4.5σ it needs
longer integrations (10–15 min/scan) and fresher baselines to be conclusive.
A fifth scan died when the SDR was unplugged mid-integration (laptop died); the
wall-clock timer kept "integrating" for 19.4 h — that bug is fixed (see
Robustness), the phantom rows were removed (backup:
`galactic_plane_log.csv.bak-20260717-rebuild`) and its npz moved to
`spectra/quarantine/`. The remaining four scans' CSV rows were re-detected with
the new smoothed algorithm. Self-test passing.
