#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dish Altitude (elevation) Angle Helper
======================================

The app's "Altitude (deg)" field is the elevation of the dish's BEAM above the
horizon. The easiest thing to physically measure is the tilt of the FEED BOOM
-- the arm that holds the LNA/feed -- because it points straight down the beam.
(0 deg = boom horizontal / looking at the horizon; 90 deg = boom straight up.)

This tool converts a couple of tape-measure numbers into that angle. Two ways:

  Method 1  (rise & run):  lay a level out horizontally from a point on the boom,
            measure how far it reaches (run) and how much the boom rises over
            that distance (rise).  angle = atan(rise / run).

  Method 2  (board & drop): rest a straight board/level of known LENGTH flat
            along the boom, then measure the vertical HEIGHT DIFFERENCE between
            its two ends.  angle = asin(drop / length).   <-- usually easiest

Use any unit as long as both numbers use the SAME unit (inches, cm, whatever).
"""
import math


def ask_float(prompt):
    while True:
        try:
            s = input(prompt).strip()
        except EOFError:
            raise SystemExit
        try:
            return float(s)
        except ValueError:
            print("   (enter a number, e.g. 12.5)")


def main():
    print("=" * 58)
    print("  DISH ALTITUDE ANGLE HELPER")
    print("  Measure along the feed boom (the arm to the LNA);")
    print("  it points where the dish is looking.")
    print("=" * 58)
    print("  Method 1 = rise & run")
    print("  Method 2 = board length & end-to-end height drop  (easiest)")
    try:
        m = input("Pick method [1/2] (default 2): ").strip()
    except EOFError:
        m = "2"
    if m == "":
        m = "2"

    if m == "1":
        run = ask_float("  Horizontal run (out from the boom): ")
        rise = ask_float("  Vertical rise over that run: ")
        if run <= 0:
            print("  Run must be positive."); return
        angle = math.degrees(math.atan2(rise, run))
    else:
        length = ask_float("  Board / boom length: ")
        drop = ask_float("  Height difference between the two ends: ")
        if length <= 0 or abs(drop) > length:
            print("  Check numbers: the height drop can't exceed the length.")
            return
        angle = math.degrees(math.asin(drop / length))

    print()
    print("  " + "-" * 40)
    print(f"    >>>  Dish Altitude = {angle:.1f} degrees  <<<")
    print("  " + "-" * 40)
    print("  Type that into the app's 'Altitude (deg)' field.")
    print()
    print("Reference (target galactic longitude -> altitude, Long Island):")
    for lon, alt in [(10, 29), (20, 38), (30, 47), (40, 56),
                     (50, 64), (60, 73), (70, 82), (80, 90)]:
        print(f"    l = {lon:>2} deg  ->  altitude {alt} deg")
    print()
    print("Tip: to double-check your mounting, point at the SUN when it is due")
    print("south -- it is a strong 1420 MHz source. The altitude that maximizes")
    print("the signal should match the Sun's known altitude at that moment.")


if __name__ == "__main__":
    main()
