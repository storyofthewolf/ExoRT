"""
check_load_rtprofile.py -- load_rtprofile must read the same variables as
input_profile in io_1D.F90.

load_rtprofile (exort_api.py) mirrors the Fortran reader so that a profile
gives the same column through exort.exe and through libexort. When a variable
is added to input_profile but not to the Python lists, the library silently
runs without it (as happened with cfrc, grav and scon). This compares the
names of the req_* / opt_* reads in io_1D.F90 with load_rtprofile's required
and optional lists, statically: no library, NetCDF file or cffi is needed.

Usage:

    python check_load_rtprofile.py      # exit 0 = PASS, 1 = FAIL
"""

import ast
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
IO_1D = os.path.join(HERE, "..", "..", "source", "src.main", "io_1D.F90")
API = os.path.join(HERE, "exort_api.py")


def fortran_reads():
    """{'required': set, 'optional': set} from the req_*/opt_* calls."""
    src = open(IO_1D).read()
    reads = {"required": set(), "optional": set()}
    for kind, name in re.findall(r"\b(req|opt)_[a-z_]+\('(\w+)'", src):
        reads["required" if kind == "req" else "optional"].add(name)
    return reads


def python_reads():
    """{'required': set, 'optional': set} from load_rtprofile's lists; the
    names it reads directly with read("...") (mw, cp) count as required."""
    tree = ast.parse(open(API).read())
    func = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "load_rtprofile")
    reads = {"required": set(), "optional": set()}
    for node in ast.walk(func):
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id in reads):
            reads[node.targets[0].id] |= {e.value for e in node.value.elts}
        elif (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
              and node.func.id == "read" and isinstance(node.args[0], ast.Constant)):
            reads["required"].add(node.args[0].value)
    return reads


def main():
    f, p = fortran_reads(), python_reads()
    ok = True
    for kind in ("required", "optional"):
        missing, extra = f[kind] - p[kind], p[kind] - f[kind]
        print(f"{kind:9s}: io_1D.F90 {len(f[kind])}, load_rtprofile {len(p[kind])}")
        if missing:
            print(f"  read by io_1D.F90 but not by load_rtprofile: {sorted(missing)}")
            ok = False
        if extra:
            print(f"  read by load_rtprofile but not by io_1D.F90: {sorted(extra)}")
            ok = False
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
