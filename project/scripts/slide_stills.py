"""Extract one representative still frame from each presentation animation.

The Beamer deck (presentation/deck.tex) is a PDF, and a PDF cannot play an
MP4. The recording workflow in presentation/SHOTLIST.md already plays each
animation fullscreen as its own take and splices it in during the edit, so
the slide only needs a single frame as a visual anchor for the viewer and
for anyone reading the PDF on its own.

Frames are pulled straight from the MP4s in results/ (not the GIF twins,
which are 256-colour quantised) using the ffmpeg binary that ships inside
the imageio-ffmpeg wheel, so nothing beyond scripts/requirements.txt is
needed. Each clip names the fraction of its runtime to sample, chosen so
the frame shows the animation's point rather than its empty first second.
"""

import argparse
import pathlib
import sys

import imageio.v2 as imageio

# (source mp4 under results/, output name, fraction of the clip's runtime).
SHOTS = [
    ("demo_timevalue.mp4", "still_timevalue.png", 0.55),
    # Early in the sweep xi is still large and the smile is at its most
    # curved; by the end it has flattened into Black-Scholes.
    ("xi_smile.mp4", "still_smile.png", 0.08),
    ("bs_collapse.mp4", "still_collapse.png", 0.90),
    ("memory_cache.mp4", "still_cache.png", 0.60),
]


def grab_frame(source: pathlib.Path, fraction: float):
    """Return the frame at `fraction` of the clip, or the last one if the
    clip turns out shorter than its metadata claimed."""
    reader = imageio.get_reader(source)
    meta = reader.get_meta_data()
    total = int(meta.get("duration", 0) * meta.get("fps", 0))
    target = max(0, int(total * fraction) - 1)
    frame = None
    for index, frame in enumerate(reader):
        if index >= target:
            break
    reader.close()
    return frame


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", default="results", type=pathlib.Path)
    parser.add_argument("--out", default="presentation/stills",
                        type=pathlib.Path)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    missing = 0
    for clip, name, fraction in SHOTS:
        source = args.results / clip
        if not source.exists():
            # Not fatal: the deck still compiles, the frame just shows the
            # draft box. Rendering the animations first fixes it.
            print(f"skip: {source} not found (make demo-anim etc. first)",
                  file=sys.stderr)
            missing += 1
            continue
        frame = grab_frame(source, fraction)
        imageio.imwrite(args.out / name, frame)
        print(f"wrote {args.out / name}")
    return 1 if missing == len(SHOTS) else 0


if __name__ == "__main__":
    sys.exit(main())
