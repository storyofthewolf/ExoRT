#!/usr/bin/env python3
"""
populate3Dmodels.py — regenerate and verify the 3dmodels/ CAM bundles from source/.

Each 3dmodels/src.cam.<version>/ directory is a *derived artifact*: the flat union
of a version directory (source/src.<version>/) and the CAM-relevant subset of
source/src.main/, which is the single directory CESM/ExoCAM consumes via -usr_src.
Because the 1-D build (build/Makefile) reads only source/, the 3dmodels/ copies are
pure duplicates and drift silently. This tool makes that duplication safe:

  check          diff every bundle file against its source origin; exit 1 on drift
  regenerate     copy source -> bundle so the bundle exactly matches source
  list           show the version->bundle mapping and per-file origin
  diff           print a unified diff of every drifted file

Reproducibility model
----------------------
The legacy bundles (n28archean, n42h2o, n68h2o) need not be committed: delete them
and recreate on demand with --regenerate when an old result must be reproduced. The
active bundles (n68equiv, n84equiv) are kept committed and guarded by --check.

Two files in a bundle are intentionally bundle-local and are NEVER synced or diffed:
  sys_rootdir.F90   differs per machine (CESM root path)
  README            bundle-specific documentation

Special-case bundles excluded from all operations
-------------------------------------------------
  src.cam7.n68equiv     CESM3 linkage; does not map cleanly to source/. Preserve as-is.
  src.cam.n68equiv.haze diverges from the primary until the haze merge lands.

These are skipped by default. Operate on them only with an explicit --include-special
(which still refuses to regenerate them — they have no faithful source origin).
"""

import argparse
import difflib
import filecmp
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Layout: this script lives in <repo>/tools/, so the repo root is its parent.
# ---------------------------------------------------------------------------
REPO      = Path(__file__).resolve().parent.parent
SOURCE    = REPO / "source"
MAIN_DIR  = SOURCE / "src.main"
BUNDLES   = REPO / "3dmodels"

# ---------------------------------------------------------------------------
# The CAM bundle manifest (exact; matches build/Makefile object lists minus the
# 1-D-only entry/IO files and all of src.misc, which CESM provides).
#
#   MAIN_FILES     come from source/src.main/
#   VERSION_FILES  come from source/src.<version>/
#   LOCAL_FILES    are bundle-local: never synced, never diffed
# ---------------------------------------------------------------------------
MAIN_FILES = [
    "exo_init_ref.F90",
    "exo_radiation_mod.F90",
    "exo_radiation_cam_intr.F90",
    "mcica.F90",
    "mcica_random_numbers.F90",
    "planck_mod.F90",
    "rayleigh_data.F90",
]

VERSION_FILES = [
    "calc_opd_mod.F90",
    "cloud.F90",
    "kabs.F90",
    "initialize_rad_mod_cam.F90",
    "model_specific.F90",
    "rad_interp_mod.F90",
    "radgrid.F90",
    "spectral_output_cam.F90",
]

# Bundle-local files: present in the bundle, deliberately NOT mirrored from source.
LOCAL_FILES = ["sys_rootdir.F90", "README"]

# version key -> (source version dir name, bundle dir name)
VERSIONS = {
    "n28archean": ("src.n28archean", "src.cam.n28archean"),
    "n42h2o":     ("src.n42h2o",     "src.cam.n42h2o"),
    "n68h2o":     ("src.n68h2o",     "src.cam.n68h2o"),
    "n68equiv":   ("src.n68equiv",   "src.cam.n68equiv"),
    "n84equiv":   ("src.n84equiv",   "src.cam.n84equiv"),
}

# Bundles that have no faithful source origin. Skipped unless --include-special,
# and never regenerated even then.
SPECIAL_BUNDLES = ["src.cam7.n68equiv", "src.cam.n68equiv.haze"]

LEGACY = {"n28archean", "n42h2o", "n68h2o"}


def manifest(version):
    """Yield (origin_path, bundle_relname) pairs for one version's synced files."""
    src_dir = SOURCE / VERSIONS[version][0]
    for f in MAIN_FILES:
        yield MAIN_DIR / f, f
    for f in VERSION_FILES:
        yield src_dir / f, f


def bundle_dir(version):
    return BUNDLES / VERSIONS[version][1]


# ---------------------------------------------------------------------------
# Verbs
# ---------------------------------------------------------------------------
def do_list(versions):
    print(f"repo root: {REPO}")
    print(f"\n{'version':<12} {'source dir':<18} -> bundle")
    print("-" * 60)
    for v in versions:
        srcname, bundlename = VERSIONS[v]
        tag = "  [legacy]" if v in LEGACY else ""
        present = "" if bundle_dir(v).is_dir() else "  (bundle absent)"
        print(f"{v:<12} {srcname:<18} -> 3dmodels/{bundlename}{tag}{present}")
    print(f"\nper-bundle files: {len(MAIN_FILES)} from src.main + "
          f"{len(VERSION_FILES)} version-specific")
    print(f"bundle-local (never synced): {', '.join(LOCAL_FILES)}")
    print(f"special bundles (skipped): {', '.join(SPECIAL_BUNDLES)}")


def do_check(versions):
    """Return number of drifted/missing files across all requested versions."""
    problems = 0
    for v in versions:
        bdir = bundle_dir(v)
        print(f"\n<<< {v}  (3dmodels/{VERSIONS[v][1]}) >>>")
        if not bdir.is_dir():
            print("  bundle directory absent — run --regenerate to create it")
            problems += 1
            continue
        for origin, name in manifest(v):
            target = bdir / name
            if not origin.exists():
                print(f"  MISSING SOURCE  {name}  (expected {origin})")
                problems += 1
            elif not target.exists():
                print(f"  MISSING BUNDLE  {name}")
                problems += 1
            elif not filecmp.cmp(origin, target, shallow=False):
                print(f"  DRIFT           {name}")
                problems += 1
            else:
                print(f"  ok              {name}")
    return problems


def do_diff(versions):
    """Print a unified diff (source -> bundle) for every drifted file.

    Returns the number of drifted/missing files (same accounting as do_check),
    so it can gate an exit status, but its purpose is to show *what* changed.
    """
    problems = 0
    for v in versions:
        bdir = bundle_dir(v)
        if not bdir.is_dir():
            print(f"<<< {v}  (3dmodels/{VERSIONS[v][1]}) >>>")
            print("  bundle directory absent — run regenerate to create it\n")
            problems += 1
            continue
        for origin, name in manifest(v):
            target = bdir / name
            if not origin.exists():
                print(f"### {v}/{name}: MISSING SOURCE (expected {origin})\n")
                problems += 1
            elif not target.exists():
                print(f"### {v}/{name}: MISSING BUNDLE\n")
                problems += 1
            elif not filecmp.cmp(origin, target, shallow=False):
                problems += 1
                # Some bundle files carry stray non-UTF-8 bytes; decode
                # leniently so the diff still renders instead of crashing.
                fromtext = origin.read_text(errors="replace")
                totext   = target.read_text(errors="replace")
                fromlines = fromtext.splitlines(keepends=True)
                tolines   = totext.splitlines(keepends=True)
                diff = difflib.unified_diff(
                    fromlines, tolines,
                    fromfile=f"source: {origin.relative_to(REPO)}",
                    tofile=f"bundle: {target.relative_to(REPO)}")
                sys.stdout.writelines(diff)
                print()  # blank line between files
    return problems


def do_regenerate(versions):
    for v in versions:
        bdir = bundle_dir(v)
        bdir.mkdir(parents=True, exist_ok=True)
        print(f"\n<<< regenerate {v} -> 3dmodels/{VERSIONS[v][1]} >>>")
        for origin, name in manifest(v):
            if not origin.exists():
                print(f"  SKIP (no source)  {name}")
                continue
            (bdir / name).write_bytes(origin.read_bytes())
            print(f"  wrote             {name}")
        # Warn if bundle-local files are absent (they are not generated here).
        for name in LOCAL_FILES:
            if not (bdir / name).exists():
                print(f"  note: bundle-local {name} absent — create it by hand "
                      f"(machine path / docs)")


# ---------------------------------------------------------------------------
def add_selection_args(sub):
    """Attach the shared version-selection flags to a subparser."""
    sel = sub.add_argument_group(
        "version selection (default: all mapped versions)")
    for v in VERSIONS:
        sel.add_argument(f"--{v}", action="store_true", help=f"operate on {v}")
    sel.add_argument("--legacy", action="store_true",
                     help="select the legacy trio (n28archean, n42h2o, n68h2o)")
    sel.add_argument("--include-special", action="store_true",
                     help=f"acknowledge {', '.join(SPECIAL_BUNDLES)} "
                          "(reported, never regenerated)")


def main():
    p = argparse.ArgumentParser(
        description="Regenerate/verify 3dmodels CAM bundles from source/.")
    sub = p.add_subparsers(dest="command", required=True, metavar="command")

    add_selection_args(sub.add_parser(
        "check", help="diff bundles against source; exit 1 on any drift"))
    add_selection_args(sub.add_parser(
        "diff", help="print a unified diff of every drifted file"))
    add_selection_args(sub.add_parser(
        "regenerate", help="copy source -> bundle so the bundle matches source"))
    add_selection_args(sub.add_parser(
        "list", help="show version->bundle mapping and exit"))

    args = p.parse_args()

    chosen = [v for v in VERSIONS if getattr(args, v)]
    if args.legacy:
        chosen += [v for v in LEGACY if v not in chosen]
    if not chosen:
        chosen = list(VERSIONS)

    if args.include_special:
        print("note: special bundles have no faithful source origin and are "
              "reported only:")
        for s in SPECIAL_BUNDLES:
            present = (BUNDLES / s).is_dir()
            print(f"      {s}  {'(present)' if present else '(absent)'}")

    if args.command == "list":
        do_list(chosen)
        return 0
    if args.command == "regenerate":
        do_regenerate(chosen)
        return 0
    if args.command == "diff":
        n = do_diff(chosen)
        return 1 if n else 0
    if args.command == "check":
        n = do_check(chosen)
        print(f"\n{'=' * 60}")
        if n:
            print(f"FAIL: {n} file(s) drifted or missing. "
                  f"Run 'regenerate' to resync (after confirming the source "
                  f"side is the intended truth).")
            return 1
        print("PASS: all selected bundles match source.")
        return 0


if __name__ == "__main__":
    sys.exit(main())
