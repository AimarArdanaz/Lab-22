#!/usr/bin/env python3
"""
measure_gel.py — measures the bands of an agarose gel from its photo.

What it does, in order:
  1. Finds the wells (the edge of the comb marks the start of each lane).
  2. Extracts the intensity profile of each lane, subtracting the local
     background, and detects bands as peaks of that profile.
  3. If you tell it which lanes are the marker and which ladder you used, it
     fits a line log10(size) = a + b * (migration distance) and uses it to
     estimate the size of the sample bands.
  4. Optionally saves an enhanced version of the image for a figure.

Why a log-linear fit: in an agarose gel, within the size range it resolves
well, the distance a fragment migrates is roughly proportional to the
logarithm of its size. If the R² of the fit is low, the marker bands have
been assigned wrongly or the gel is compressed at the ends: do not trust the
estimated sizes.

The row of wells is found automatically (it tries horizontal strips and keeps
the one that gives the most evenly spaced wells). If it gets it wrong, pass it
by hand with --wells-y. If the photo has two gels or two rows of wells,
analyse one at a time: the one it picks by itself, and the other one by
passing its --wells-y.

It also detects whether the photo is a negative (dark bands on a light
background) and inverts it; this can be forced with --invert yes/no.

Usage:
    python measure_gel.py photo.jpg                        # detect bands only
    python measure_gel.py photo.jpg --markers 1 4 \
        --ladder 5000 4000 3000 2000 1500 1000 700 500
    python measure_gel.py photo.jpg --enhance output.png

Requirements: Python 3 with numpy and Pillow.
"""
import argparse
import itertools
import numpy as np
from PIL import Image

# Known ladders: sizes in bp from largest to smallest, and which ones are the
# reference bands (the ones the manufacturer makes brighter for orientation).
LADDERS = {
    # Thermo Scientific, SM1331/SM1332 and the ready-to-use SM1333/SM1334.
    # Manual MAN0013047 Rev. D.00: "fifteen ... individual DNA fragments (in base
    # pairs): 20000, 10000, 7000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500,
    # 400, 300, 200, 75. It contains three reference bands (5000, 1500 and 500 bp)".
    "generuler-1kb-plus": {
        "sizes": [20000, 10000, 7000, 5000, 4000, 3000, 2000, 1500,
                  1000, 700, 500, 400, 300, 200, 75],
        "reference": [5000, 1500, 500],
        "name": "GeneRuler 1 kb Plus (Thermo SM1331/SM1333)",
    },
    "generuler-1kb": {
        "sizes": [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500,
                  2000, 1500, 1000, 750, 500, 250],
        "reference": [6000, 3000, 1000],
        "name": "GeneRuler 1 kb (Thermo SM0311/SM0313)",
    },
}
DEFAULT_LADDER = "generuler-1kb-plus"


def load(path):
    return np.asarray(Image.open(path).convert("L"), dtype=float)


def is_negative(a):
    """Are the bands light on a dark background, or the other way round?

    Looks only at the central area (the gel, without the black frame of the
    photo) and compares how far a band sticks out above and below the median
    of the background. Returns True if the image has to be inverted.
    """
    h, w = a.shape
    c = a[int(h * .2):int(h * .8), int(w * .2):int(w * .8)]
    med = np.median(c)
    return (med - np.percentile(c, 0.2)) > (np.percentile(c, 99.8) - med)


def find_wells(a, y0=150, y1=210, min_width=18):
    """Wells are rectangles with an edge: they show up as a peak of vertical
    gradient. Returns the x centre of each one."""
    strip = a[y0:y1, :]
    g = np.abs(np.diff(strip, axis=0)).sum(axis=0)
    g = np.convolve(g, np.ones(5) / 5, mode="same")
    threshold = np.percentile(g, 60)
    centres, inside, start = [], False, 0
    for x, v in enumerate(g):
        if v > threshold and not inside:
            start, inside = x, True
        elif v <= threshold and inside:
            if x - start > min_width:
                centres.append((start + x) / 2)
            inside = False
    return centres


def find_wells_auto(a, strip_height=55, step=8, min_width=18):
    """Tries horizontal strips across the whole upper half of the photo and
    keeps the one that gives the most wells, best aligned.

    The wells of a comb are evenly spaced, so the right strip is the one that
    gives many centres with almost identical spacings.
    Score = number of wells / (1 + relative variation of the spacings).
    Returns (centres, y0, y1, score).
    """
    height = a.shape[0]
    best = ([], 0, strip_height, 0.0)
    for y0 in range(20, int(height * 0.55), step):
        y1 = y0 + strip_height
        if y1 >= height:
            break
        c = find_wells(a, y0, y1, min_width)
        if len(c) < 3:
            continue
        sep = np.diff(c)
        cv = sep.std() / sep.mean() if sep.mean() else 9.9
        score = len(c) / (1 + 5 * cv)
        if score > best[3]:
            best = (c, y0, y1, score)
    return best


def lane_profile(a, xc, width=46, y0=0, y1=None, background_window=40):
    """Vertical profile of a lane, with the local background subtracted."""
    y1 = y1 or a.shape[0]
    col = a[y0:y1, int(xc - width / 2):int(xc + width / 2)].mean(axis=1)
    background = np.array([np.median(col[max(0, i - background_window):i + background_window])
                           for i in range(len(col))])
    return col - background


def find_peaks(p, min_height=4.0, separation=7):
    """Local maxima of the profile: each one is a band.

    A wide band can give two maxima side by side (a plateau), so at the end
    peaks closer than `separation` are merged and the highest one is kept.
    Without this, the same band would be counted twice.
    """
    raw = []
    for y in range(len(p)):
        v = p[y]
        if v >= min_height and v == p[max(0, y - separation):y + separation + 1].max():
            raw.append((y, float(v)))
    out = []
    for y, v in raw:
        if out and y - out[-1][0] < separation:
            if v > out[-1][1]:
                out[-1] = (y, v)
        else:
            out.append((y, v))
    return out


def fit_loglinear(distances, sizes):
    """Line log10(size) versus migration distance. Returns (a, b, R²)."""
    d = np.asarray(distances, float)
    y = np.log10(np.asarray(sizes, float))
    b, a = np.polyfit(d, y, 1)
    pred = a + b * d
    ss_res = ((y - pred) ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    return a, b, 1 - ss_res / ss_tot


def assign_ladder(marker_bands, ladder, reference, max_comb=400000):
    """Decides which marker band corresponds to each size of the ladder.

    The problem: a gel almost never resolves all 15 bands of a GeneRuler.
    Some are lost through exposure and others merge at the top, where the gel
    does not separate. So pairing "the first band I see with the largest size"
    does not work.

    The clue that solves it is brightness. The manufacturer deliberately loads
    some reference bands with more DNA for orientation (in the GeneRuler 1 kb
    Plus, 5000, 1500 and 500 bp). So **the three brightest bands of the marker
    lane must be those three**, in that order from top to bottom. That fixes
    three points of the line and leaves little room for the rest.

    With that constraint, every way of distributing the remaining bands among
    the remaining sizes (keeping the order) is tried, and the one that best
    fits the line log10(size) = a + b·distance wins.

    If nothing comes out with that constraint, it tries without it and warns:
    that usually means the lane is not a marker, reference bands are missing
    through exposure, or the declared ladder is not the one that was loaded.

    Returns (list of (y, height, size), R², number of anchored references,
    mismatch) or None.
    """
    n = len(marker_bands)
    if n < 3 or n > len(ladder):
        return None
    d = np.array([y for y, _ in marker_bands], float)
    h = np.array([hh for _, hh in marker_bands], float)

    # positions of the brightest peaks, sorted from top to bottom
    n_ref = min(len(reference), n)
    anchors = {}
    if reference:
        brightest = sorted(sorted(range(n), key=lambda k: -h[k])[:n_ref])
        anchors = {k: reference[j] for j, k in enumerate(brightest)}

    def search(anchors):
        best, tried = None, 0
        for idx in itertools.combinations(range(len(ladder)), n):
            tried += 1
            if tried > max_comb:
                break
            sizes = [ladder[i] for i in idx]
            if any(sizes[k] != v for k, v in anchors.items()):
                continue
            y = np.log10(np.array(sizes, float))
            b, a0 = np.polyfit(d, y, 1)
            if b >= 0:                  # must go down: larger fragments migrate less
                continue
            pred = a0 + b * d
            ss_res = ((y - pred) ** 2).sum()
            ss_tot = ((y - y.mean()) ** 2).sum()
            r2 = 1 - ss_res / ss_tot if ss_tot else 0
            if best is None or r2 > best[0]:
                best = (r2, sizes)
        return best

    res = search(anchors)
    anchored = len(anchors)
    if res is None and anchors:
        res = search({})
        anchored = 0
    if res is None:
        return None
    r2, sizes = res
    mismatch = 0 if anchored == len(reference) else len(reference) - anchored
    chosen = [(int(d[k]), h[k], int(sizes[k])) for k in range(n)]
    return chosen, r2, anchored, mismatch


def enhance(a, crop, scale, gamma, dest):
    y0, y1, x0, x1 = crop
    c = a[y0:y1, x0:x1]
    lo, hi = np.percentile(c, 1), np.percentile(c, 99.8)
    c = np.clip((c - lo) / (hi - lo), 0, 1) ** gamma
    im = Image.fromarray((c * 255).astype("uint8"))
    im = im.resize((int((x1 - x0) * scale), int((y1 - y0) * scale)), Image.LANCZOS)
    im.save(dest)
    return dest


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("image")
    ap.add_argument("--wells-y", nargs=2, type=int, default=None,
                    metavar=("Y0", "Y1"),
                    help="strip in which to look for the wells; found automatically by default")
    ap.add_argument("--lanes", nargs="*", type=float,
                    help="x centres of the lanes (otherwise detected automatically)")
    ap.add_argument("--y-range", nargs=2, type=int, metavar=("Y0", "Y1"),
                    help="only look between these two heights of the photo. Useful to "
                         "leave out the dye front and the compressed top part, which "
                         "does not resolve")
    ap.add_argument("--min-height", type=float, default=4.0,
                    help="minimum height of a peak to count it as a band; lower it "
                         "if the gel is underexposed and weak bands do not show up")
    ap.add_argument("--markers", nargs="*", type=int, default=[],
                    help="lane numbers (1 = the first) that carry the marker")
    ap.add_argument("--ladder", nargs="*", default=[DEFAULT_LADDER],
                    help="name of a known ladder (" + ", ".join(LADDERS) +
                         ") or the sizes in bp by hand, from largest to smallest. "
                         f"Default: {DEFAULT_LADDER}")
    ap.add_argument("--assign", choices=["auto", "order"], default="auto",
                    help="auto: works out which size each detected band is (the usual "
                         "case, because the gel does not resolve all 15 bands). order: "
                         "pairs them one by one in the given order, only if you know "
                         "they are all visible")
    ap.add_argument("--invert", choices=["auto", "yes", "no"], default="auto",
                    help="for negative photos (dark bands on a light background)")
    ap.add_argument("--enhance", metavar="OUTPUT.png",
                    help="save an enhanced version for a figure")
    ap.add_argument("--crop", nargs=4, type=int, metavar=("Y0", "Y1", "X0", "X1"))
    ap.add_argument("--scale", type=float, default=3.2)
    ap.add_argument("--gamma", type=float, default=0.6)
    args = ap.parse_args()

    a = load(args.image)
    height, width = a.shape
    print(f"image {width}x{height}")

    invert = is_negative(a) if args.invert == "auto" else (args.invert == "yes")
    if invert:
        a = 255.0 - a
        print("negative photo: inverted so that bands are peaks")

    if args.wells_y:
        strip = tuple(args.wells_y)
        detected_lanes = find_wells(a, *strip)
    else:
        detected_lanes, sy0, sy1, score = find_wells_auto(a)
        strip = (sy0, sy1)
        print(f"well strip chosen automatically: y = {sy0}–{sy1}  (score {score:.1f})")
        if score < 2:
            print("    [warning] the row of wells is not clear. Look at the photo and "
                  "pass --wells-y Y0 Y1 by hand")

    lanes = args.lanes or detected_lanes
    if not lanes:
        raise SystemExit("no wells found: pass --lanes X1 X2 ... by hand")
    sep = np.diff(lanes)
    print(f"{len(lanes)} wells detected at x = " +
          ", ".join(f"{c:.0f}" for c in lanes))
    if len(sep) and sep.std() / sep.mean() > 0.15:
        print("    [warning] the spacing between wells is irregular: some may be "
              "missing or extra. Check the photo")
    y_well = sum(strip) / 2

    y_start = int(y_well) + 20
    y_end = height
    if args.y_range:
        y_start, y_end = max(y_start, args.y_range[0]), args.y_range[1]
        print(f"only bands between y = {y_start} and y = {y_end} are considered")

    bands = {}
    for i, xc in enumerate(lanes, 1):
        p = lane_profile(a, xc, y0=y_start, y1=y_end)
        b = find_peaks(p, args.min_height)
        bands[i] = [(y + y_start, hh) for y, hh in b]
        if b:
            print(f"\nlane {i} (x≈{xc:.0f})")
            for y, hh in bands[i]:
                near = "   [right next to the well: may be the edge, not a band]" \
                       if y - strip[1] < 15 else ""
                print(f"    y = {y:4d}   height = {hh:6.2f}   "
                      f"migration = {y - y_well:5.1f} px{near}")

    if args.markers:
        # --ladder accepts a name from the catalogue or the sizes by hand
        if len(args.ladder) == 1 and not args.ladder[0].replace(".", "").isdigit():
            key = args.ladder[0].lower()
            if key not in LADDERS:
                raise SystemExit(f"unknown ladder: {key}. "
                                 f"Known: {', '.join(LADDERS)}")
            lad = LADDERS[key]
            sizes, reference = lad["sizes"], lad["reference"]
            print(f"\nladder: {lad['name']}")
            print(f"    {len(sizes)} bands, reference (brightest) at "
                  + ", ".join(f"{r} bp" for r in reference))
        else:
            sizes = [float(x) for x in args.ladder]
            reference = []
            print(f"\nladder given by hand: {len(sizes)} bands")

        d, t_bp, results = [], [], []
        for n in args.markers:
            bs = bands.get(n, [])
            if len(bs) < 3:
                print(f"[warning] lane {n} only has {len(bs)} bands: not used. "
                      f"Lower --min-height to see the weak ones")
                continue
            if args.assign == "order":
                if len(bs) != len(sizes):
                    print(f"[warning] lane {n} has {len(bs)} bands and the ladder "
                          f"declares {len(sizes)}: it cannot be used with --assign order")
                    continue
                chosen = [(y, hh, s) for (y, hh), s in zip(bs, sizes)]
                r2_lane, anchored, mismatch = None, None, 0
            else:
                res = assign_ladder(bs, sizes, reference)
                if res is None:
                    print(f"[warning] lane {n} has {len(bs)} bands and they cannot be "
                          f"distributed among the {len(sizes)} sizes of the ladder")
                    continue
                chosen, r2_lane, anchored, mismatch = res

            print(f"\nlane {n} — assignment" +
                  (f" (R² = {r2_lane:.4f}, {anchored} of {len(reference)} reference "
                   f"bands anchored by brightness)" if r2_lane is not None else ""))
            for y, hh, s in chosen:
                mark = "  <- reference" if s in reference else ""
                print(f"    y = {int(y):4d}   height = {hh:6.2f}   ->  {int(s):6d} bp{mark}")
            if r2_lane is not None and reference and anchored < len(reference):
                print(f"    [warning] the reference bands could not be anchored by "
                      f"brightness: this assignment is a guess. Look at the photo")
            reliable = (not reference) or (mismatch == 0) or r2_lane is None
            results.append((n, chosen, r2_lane, anchored, reliable))

        good = [r for r in results if r[4]]
        if reference and good and len(good) < len(results):
            discarded = [str(r[0]) for r in results if not r[4]]
            print(f"\n[warning] only lane "
                  f"{', '.join(str(r[0]) for r in good)} is used for calibration: in lane "
                  f"{', '.join(discarded)} the bright bands do not match the "
                  f"reference ones, "
                  f"so its assignment cannot be trusted")
            results = good
        for _n, chosen, _r2, _an, _ok in results:
            for y, _hh, s in chosen:
                d.append(y - y_well); t_bp.append(s)

        if d:
            a0, b0, r2 = fit_loglinear(d, t_bp)
            print(f"\nlog-linear fit with {len(d)} marker bands: R² = {r2:.4f}")
            if r2 < 0.98:
                print("    [warning] low R²: do not trust the sizes below")
            print("\nestimated sizes")
            tmin, tmax = min(t_bp), max(t_bp)
            for i in bands:
                if i in args.markers:
                    continue
                for y, hh in bands[i]:
                    size = 10 ** (a0 + b0 * (y - y_well))
                    outside = ""
                    if size > tmax or size < tmin:
                        outside = "   [outside the calibrated range: weak estimate]"
                    print(f"    lane {i:2d}   y = {int(y):4d}   ≈ {size:8.0f} bp{outside}")

    if args.enhance:
        crop = args.crop or [0, height, 0, width]
        print("\nenhanced ->", enhance(a, crop, args.scale, args.gamma, args.enhance))


if __name__ == "__main__":
    main()
