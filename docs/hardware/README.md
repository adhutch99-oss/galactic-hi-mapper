# Photographs

| File | What it shows | Honest description |
|---|---|---|
| `bench-setup.jpg` | The whole signal chain assembled indoors, laptop running the app | Bench setup — **not** an observing session |
| `app-live-hardware.jpg` | The app streaming from the SDR, backend on Hardware, maps reloaded from the log | Off-plane at the moment of the photo (banner red) |
| `app-plane-banner.jpg` | The pointing banner close up, amber at b = +2.87° | The app in **Simulation** mode |

Still to come: the dish outdoors on its mount, and a mid-scan shot with the
banner green next to the log row it produced.

**Before committing a new photo, strip its EXIF.** Phone photos embed GPS
coordinates accurate to a few metres — the originals here did. Re-save the
pixels through an image editor; cropping in a phone gallery does not remove it.

```bash
python -c "from PIL import Image; e=Image.open('new.jpg').getexif(); print(len(e), bool(e.get_ifd(0x8825)))"
```

Both numbers must be `0` and `False`.

## Still missing

1. **The dish outdoors on its mount**, feed boom visible. The only shot that shows the
   instrument deployed rather than assembled on a floor.
2. **The app mid-scan with the banner green** and a line in the spectrum, paired with the
   log row it produced. That combination is evidence rather than illustration.
3. **The feed and LNA close up**, SAWbird power LED lit.
4. **Elevation being measured** — phone inclinometer on the feed boom.

Lowercase hyphenated filenames, under ~2 MB each.
