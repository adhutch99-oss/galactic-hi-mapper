# Graphs Guide

How to take the data, how to turn the log into plots, and how to read each one.
Companion to [OPERATING_GUIDE.md](OPERATING_GUIDE.md). Kept as fixed-width text
because the diagrams and tables rely on alignment.

```text
================================================================================
 GRAPHS GUIDE — how to make the plots, how to read them, what to do with them
================================================================================
Companion to OPERATING_GUIDE.md. This covers the whole run: taking the data,
turning the log into graphs, and turning the graphs into science.


--------------------------------------------------------------------------------
 PART 1 — TAKING THE DATA (a full observing session, step by step)
--------------------------------------------------------------------------------

 A. FIRST, KNOW YOU'RE ON REAL HARDWARE (once the parts arrive)

    Plug the SDR into USB BEFORE launching the app, then check these four
    signs, in order:

    1. The SAWbird LNA's power LED lights up as soon as the SDR is plugged
       in (the SDR's always-on bias-tee is feeding it). No LED = check the
       cable between LNA and SDR.
    2. Launch the app. The status bar (bottom left) must say:
           "RTL-SDR detected. Hardware backend available."
       and the Backend dropdown will be set to "Hardware".
       If instead it says "No SDR found — falling back to Simulation mode",
       the software did NOT find the SDR — usually the Zadig/WinUSB driver
       step was missed (OPERATING_GUIDE.md section 4).
    3. Press START LIVE PREVIEW. The status bar must say:
           "Streaming from Hardware"        (not "from Simulation")
    4. Every row you record will show [Hardware] in the status message, and
       the CSV's Source column will say Hardware. Simulation rows say
       Simulation — they can never be confused.

    Reality check: block the dish (point at the ground) — the noise floor
    should change. A simulation never reacts to the real world.

 B. ONCE PER SESSION — the cold-sky baseline (NOT a recorded measurement)

    You do NOT need to check or uncheck ANY box for this. The whole thing:

    1. Point the dish somewhere high and empty, AWAY from the Milky Way band.
       NEVER at the sun. EASIEST METHOD: keep the dish due south (don't
       rotate it), crank the altitude HIGH, and read the banner — when it
       shows |b| bigger than ~30 degrees (deep red), that IS verified cold
       sky. The app can only compute b for a south-facing dish, so staying
       south lets it confirm your baseline pointing for you. Keep trees and
       buildings out of the beam (warm ground raises the baseline).
    2. START LIVE PREVIEW, let the line settle ~10 seconds.
    3. Click "Capture Cold Sky Baseline" — ONE click. The app ticks
       "Subtract baseline" for you automatically.
    4. Leave "Subtract baseline" CHECKED for the rest of the night.
       Nothing gets unchecked, ever, during a session.

    WHAT THE GRAPH LOOKS LIKE AFTER CAPTURE: the curved bandpass disappears
    and you get a nearly FLAT noisy line sitting near 0 dB. That is correct —
    flat = "no difference from cold sky". When you aim at the galactic plane,
    the hydrogen line pokes UP out of that flat line. If the plot ever looks
    empty, the y-axis just needs a rescale — toggling any checkbox re-scales
    it (and the app now does this automatically when the view changes).

    Do NOT press LOG DRIFT SCAN or STOP & RECORD for the baseline.

 C. PLAN THE NIGHT

    Click "What can I see tonight?" — it lists every longitude transiting in
    the next 12 h with the altitude to set and the clock time. Pick a handful
    5-10 degrees apart and note their times.

 D. FOR EACH TARGET (repeat this loop all night)

    1. A few minutes before its listed time, set the dish to that target's
       altitude (type the l in the planner and click "Apply altitude").
       Dish stays pointed due south the entire night, every night.
    2. Watch the banner. Red/amber = wait. GREEN (b near 0) = the plane is
       in the beam: click LOG DRIFT SCAN.
    3. Let it integrate 1-3 minutes (the counter shows elapsed seconds).
       Don't touch the dish. Longer = cleaner for faint targets (l > 60,
       outer galaxy).
    4. Click STOP & RECORD ROW. Done — spectrum saved, CSV row(s) written,
       maps updated. One scan can log several components (one per arm).
    5. Move to the next target's altitude and wait for its green banner.

 E. ORDER DOES NOT MATTER

    You can observe any longitude, any altitude, any night, in any order —
    skip around freely. Every recorded row is stamped with its own l, b,
    velocity, and time, and is completely independent of every other row.
    The maps and all four analysis plots are rebuilt from the WHOLE log every
    time, and they bin/sort by longitude themselves — the CSV's row order is
    irrelevant. Opening the CSV in Excel is equally safe in any order.

    The only "schedule" is the sky's: each longitude is only due south at its
    transit time, so WITHIN one night the sky decides what is available when
    (that's what the tonight-list is for). Gaps you leave tonight can be
    filled next week — same longitude = same altitude forever, just ~4 min
    earlier per night.


--------------------------------------------------------------------------------
 PART 2 — HOW TO MAKE THE GRAPHS (step by step)
--------------------------------------------------------------------------------

 1. Record data with the app (see OPERATING_GUIDE.md section 8). Every time you
    press STOP & RECORD ROW, the app:
      - writes one CSV row per detected velocity component, and
      - saves the scan's full spectrum to the spectra\ folder automatically.
    You don't have to do anything extra during observing.

 2. When you want plots, double-click ONE of these — no typing needed:

       Make Plots (Analyze Log).bat      uses EVERY row (sim + hardware).
                                         Use this while practicing in
                                         Simulation mode.

       Make Plots (Hardware Only).bat    ignores all Simulation rows.
                                         Use this once the real telescope
                                         is running — it is the same thing
                                         as the --hardware-only flag below.

    Either way, four PNG files appear next to the CSV:

       lv_diagram.png       longitude-velocity peak scatter + tangent envelope
       lv_diagram_2d.png    longitude-velocity INTENSITY image (full profiles)
       rotation_curve.png   rotation speed vs distance from galactic center
       faceon_map.png       top-down map of the Milky Way (the arm map)

    Re-running overwrites them with the latest data — safe to run any time.

 3. Useful command-line versions (run in the radio folder):

       python313 analyze_log.py                    normal run
       python313 analyze_log.py --show             also pop plots on screen
       python313 analyze_log.py --hardware-only    IGNORE simulation rows
       python313 analyze_log.py --min-snr 8        stricter detection cut
       python313 analyze_log.py --fits             also export a FITS table

    (python313 = C:\Users\andre\AppData\Local\Programs\Python\Python313\python.exe)

    You only need the command line for the extra flags (--min-snr, --fits);
    the two .bat files cover the everyday cases. IMPORTANT once hardware
    arrives: use the Hardware Only bat (or --hardware-only) so any leftover
    simulation practice rows never mix into your real science plots.

 4. Quality control happens automatically: rows with SNR below 5 sigma (scans
    where nothing was really detected) are dropped before plotting. If a weird
    point survives, open galactic_plane_log.csv in Excel, find the row, and
    check its SNR_sigma and Line_FWHM_kHz columns — real HI is broad (tens of
    kHz); a super-narrow "line" is RFI.

 5. Excel alternative: open the CSV, make an XY scatter of
    Galactic_Longitude_L (x) vs Doppler_Velocity_LSR_km_s (y).
    That's a hand-made l-v diagram.


--------------------------------------------------------------------------------
 PART 3 — WHAT EACH GRAPH MEANS AND WHAT TO DO WITH IT
--------------------------------------------------------------------------------

 == lv_diagram.png (longitude-velocity scatter) ==

 WHAT IT IS: every detected line peak as one dot; x = where you looked along
 the plane, y = gas velocity (LSR), color = signal strength. The dashed gray
 curve is the "tangent envelope" — the maximum velocity a flat-rotating galaxy
 allows at each longitude.

 HOW TO READ IT:
   - Inner galaxy (l = 10-80): dots should be POSITIVE velocity (redshifted).
   - A dot ON the envelope = gas at the tangent point (closest to the galactic
     center your sightline gets). These are the rotation-curve gold.
   - Several dots stacked at ONE longitude = several arms on that sightline.
   - Dots forming a diagonal LANE across many longitudes = one spiral arm.
   - Anticenter (l ~ 180): velocities near 0. Outer arc (l = 200-260): mild.

 WHAT TO DO WITH IT:
   - Look for lanes. 2-3 distinct lanes across l = 20-70 means you have
     resolved separate arms (Sagittarius-Carina, Scutum-Centaurus...).
   - Any dot far ABOVE the envelope is suspect — check its CSV row for RFI.
   - Gaps in longitude = observe those transits on the next clear night.

 == lv_diagram_2d.png (the intensity image — your best plot) ==

 WHAT IT IS: the same axes, but built from the FULL saved spectra: every
 velocity channel of every scan, not just the peaks. This is the same kind of
 figure professional HI surveys publish.

 HOW TO READ IT: bright ridges = hydrogen. A continuous bright ridge snaking
 across longitudes IS a spiral arm — including blended shoulders the
 peak-finder can't split. The "x" marks show where the logged peaks landed;
 they should sit on the ridges.

 WHAT TO DO WITH IT: this is your headline figure. Compare its ridge pattern
 with a published Milky Way l-v diagram (search: "LAB survey longitude
 velocity diagram") — matching ridge shapes is the proof your telescope is
 seeing real galactic structure.

 == rotation_curve.png (rotation speed vs radius) ==

 WHAT IT IS: from each inner-galaxy longitude, the most extreme velocity is
 assumed to come from the tangent point at radius R = 8.2 * sin(l) kpc, giving
 one measured rotation speed V(R) per longitude.

 HOW TO READ IT: x = distance from the galactic center (kpc), y = how fast gas
 orbits there (km/s). The dashed line is the textbook 220 km/s.

 WHAT TO DO WITH IT:
   - A planet-style system would show V dropping like 1/sqrt(R) at large R.
   - The real Milky Way stays roughly FLAT (~200-230 km/s) out to the Sun.
   - IMPORTANT: this plot is only a measurement if the TERMINAL velocity at
     each longitude was actually measured. The tangent-point method assumes
     the most extreme velocity in the profile comes from the tangent point.
     If the detector logged the bright local-gas peak near v = 0 instead,
     the plotted curve is just the geometric V0*sin(l) term rearranged, and
     it will look convincingly flat while measuring nothing. Check each
     profile by eye before trusting the curve.
   - Needs points spread over l = 15-75. More longitudes = better curve.
     Below ~5 longitudes it will look ragged; don't over-interpret it early.

 == faceon_map.png (top-down Milky Way map) ==

 WHAT IT IS: each (longitude, velocity) is converted to a distance (kinematic
 distance, assuming flat rotation) and drawn on a bird's-eye view. Yellow star
 = Galactic Center, blue dot = Sun, dashed circle = Sun's orbit.

 HOW TO READ IT:
   - Filled circles = "near" distance solution, open circles = "far" solution.
     Inside the Sun's orbit, one velocity honestly matches TWO distances —
     both are drawn; only one of each pair is physically real.
   - Squares = tangent points: NO ambiguity, these positions are solid.
   - Diamonds = outer galaxy: also unambiguous.

 WHAT TO DO WITH IT:
   - Trust the squares and diamonds first; arcs they trace = arm segments.
   - When near/far pairs are drawn, the one that lines up with an arc of
     tangent points is usually the real one.
   - This map IS "mapping the spiral arms of the galaxy." Overlay it mentally
     (or in an image editor) on a Milky Way artist's map — your points should
     land along the drawn arms.

 == The three LIVE graphs in the app (quick reference) ==

   Top (spectrum):   the measurement. Peak's distance LEFT of the green line
                     = redshift, read km/s directly off the TOP axis.
   Middle (l-b map): coverage bookkeeping — which sky you've already done.
   Bottom (live l-v): same as lv_diagram.png, growing in real time.


--------------------------------------------------------------------------------
 PART 4 — THE GREEN LINE vs THE GREEN BANNER (common confusion, read this!)
--------------------------------------------------------------------------------

 DO NOT wait for the spike to line up with the green dashed line. That is not
 what the green line is for, and on target it will almost never happen.

   - The GREEN BANNER (b near 0) tells you WHEN to record: the galactic plane
     is crossing your beam. That is the thing you wait for.
   - The GREEN LINE on the spectrum is the ZERO-VELOCITY reference
     (1420.4058 MHz). The GAP between your spike and the green line IS the
     Doppler measurement — the entire point of the project.

 For inner-galaxy targets (l = 10-80) a correct, on-target detection is a
 spike sitting clearly LEFT of the green line (redshifted, +30 to +120 km/s).
 A spike sitting exactly ON the green line would mean gas with zero radial
 velocity — boring, and usually just local foreground hydrogen.

 So the workflow is:  green BANNER -> start LOG DRIFT SCAN -> integrate 1-3
 min -> STOP & RECORD. The spike stays offset from the green line the whole
 time, and that offset is your data.

--------------------------------------------------------------------------------
 Checklist for a "publishable" (resume-ready) result
--------------------------------------------------------------------------------
 [ ] 8+ longitudes between l = 15 and 80 (rotation curve + inner arms)
 [ ] 4+ longitudes between l = 180 and 250 (outer/Perseus arm)
 [ ] all plots regenerated with --hardware-only
 [ ] lv_diagram_2d ridges compared against a published l-v survey figure
 [ ] terminal velocity verified by eye on each profile before the
     rotation curve is quoted as a measurement
 [ ] faceon_map arcs identified with named arms (Sagittarius, Perseus...)
================================================================================
```
