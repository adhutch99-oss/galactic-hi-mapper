# Photographs

| File | What it shows | Honest description |
|---|---|---|
| `bench-setup.jpg` | The whole signal chain assembled indoors, laptop running the app | Bench setup — **not** an observing session |
| `app-live-hardware.jpg` | The app streaming from the SDR, backend on Hardware, maps reloaded from the log | Off-plane at the moment of the photo (banner red) |
| `app-plane-banner.jpg` | The pointing banner close up, amber at b = +2.87° | The app in **Simulation** mode |

## Two rules these follow

**Nothing implies an observation that did not happen.** The bench photo is indoors with the
dish aimed at a wall. The interface photos have a red and an amber banner, and one of them
is the simulator. All three are captioned as such in the top-level README. They are worth
publishing because they show the receiver chain is real hardware and that the software
drives it — that is a different claim from "here is a detection", and the spectra in
`spectra/` carry that one.

**No location data, in the pixels or in the file.** Phone photos embed GPS coordinates in
EXIF, typically accurate to a few metres; the originals here did. Every image in this folder
was re-encoded from raw pixels with no metadata attached, and the latitude/longitude fields
visible on screen were blanked. Check any new photo before committing it:

```bash
python -c "from PIL import Image; e=Image.open('new.jpg').getexif(); print(len(e), bool(e.get_ifd(0x8825)))"
```

Both numbers must be `0` and `False`. Stripping EXIF means opening the image and re-saving
the pixels alone — renaming the file or cropping it in a phone gallery does not do it.

## Still missing

1. **The dish outdoors on its mount**, feed boom visible. The only shot that shows the
   instrument deployed rather than assembled on a floor.
2. **The app mid-scan with the banner green** and a line in the spectrum, paired with the
   log row it produced. That combination is evidence rather than illustration.
3. **The feed and LNA close up**, SAWbird power LED lit.
4. **Elevation being measured** — phone inclinometer on the feed boom.

Lowercase hyphenated filenames, under ~2 MB each.
