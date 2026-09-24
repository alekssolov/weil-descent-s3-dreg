#!/usr/bin/env python3
"""
m4gb_gm_patch.py -- fix the Gebauer-Moeller F-criterion in M4GB's update()
(src/m4gb.hpp, 2018-01-31 revision, lines ~1634-1644).

As written, a new critical pair A is discarded whenever ANY other new pair B
satisfies  lcm(B) | lcm(A)  with NON-STRICT divisibility.  Two new pairs with
EQUAL lcm therefore discard each other, whereas Gebauer-Moeller keeps exactly
one of them (none if a coprime pair shares that lcm).  The missing
S-polynomials are recovered later through pairs of higher lcm degree, which
inflates the reported max lcmdegree (e.g. 4 instead of 3 for the HKY n=16
Weil-descent systems).

Usage:  python3 m4gb_gm_patch.py /cygdrive/c/m4gb/src/m4gb.hpp
A backup  m4gb.hpp.orig  is written next to the file.
"""
import sys

OLD = (
    "\t\t\tfor (auto it = newCP_good.begin(); it != newCP_good.end(); ++it)\n"
    "\t\t\t{\n"
    "\t\t\t\tif (!it->second)\n"
    "\t\t\t\t\tcontinue;\n"
    "\t\t\t\tfor (auto it2 = newCP_good.begin(); it2 != newCP_good.end() && it2->first.lcmdeg <= it->first.lcmdeg; ++it2)\n"
    "\t\t\t\t\tif (it != it2 && it2->first.intlcm | it->first.intlcm)\n"
    "\t\t\t\t\t{\n"
    "\t\t\t\t\t\tit->second = false;\n"
    "\t\t\t\t\t\tbreak;\n"
    "\t\t\t\t\t}\n"
    "\t\t\t}\n"
)

NEW = (
    "\t\t\t// Gebauer-Moeller, corrected:\n"
    "\t\t\t//  (M) drop a pair if another new pair has an lcm that PROPERLY divides its lcm;\n"
    "\t\t\t//  (F) among new pairs with EQUAL lcm keep exactly one -- none if one of them is coprime.\n"
    "\t\t\t// (original code used non-strict divisibility, so equal-lcm pairs discarded each other)\n"
    "\t\t\t{\n"
    "\t\t\t\tstd::set<int_monomial_t> lcm_coprime, lcm_kept;\n"
    "\t\t\t\tfor (auto it = newCP_good.begin(); it != newCP_good.end(); ++it)\n"
    "\t\t\t\t\tif (!it->second)\n"
    "\t\t\t\t\t\tlcm_coprime.insert(it->first.intlcm);\n"
    "\t\t\t\tfor (auto it = newCP_good.begin(); it != newCP_good.end(); ++it)\n"
    "\t\t\t\t{\n"
    "\t\t\t\t\tif (!it->second)\n"
    "\t\t\t\t\t\tcontinue;\n"
    "\t\t\t\t\tif (lcm_coprime.count(it->first.intlcm) != 0)\n"
    "\t\t\t\t\t{\n"
    "\t\t\t\t\t\tit->second = false;\n"
    "\t\t\t\t\t\tcontinue;\n"
    "\t\t\t\t\t}\n"
    "\t\t\t\t\tfor (auto it2 = newCP_good.begin(); it2 != newCP_good.end() && it2->first.lcmdeg <= it->first.lcmdeg; ++it2)\n"
    "\t\t\t\t\t\tif (it != it2 && it2->first.intlcm != it->first.intlcm && it2->first.intlcm | it->first.intlcm)\n"
    "\t\t\t\t\t\t{\n"
    "\t\t\t\t\t\t\tit->second = false;\n"
    "\t\t\t\t\t\t\tbreak;\n"
    "\t\t\t\t\t\t}\n"
    "\t\t\t\t\tif (!it->second)\n"
    "\t\t\t\t\t\tcontinue;\n"
    "\t\t\t\t\tif (!lcm_kept.insert(it->first.intlcm).second)\n"
    "\t\t\t\t\t\tit->second = false; // an earlier pair with the same lcm is already kept\n"
    "\t\t\t\t}\n"
    "\t\t\t}\n"
)

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "src/m4gb.hpp"
    with open(path, newline="") as f:
        src = f.read()
    crlf = "\r\n" in src
    text = src.replace("\r\n", "\n")
    if NEW in text:
        print("already patched:", path)
        sys.exit(0)
    if text.count(OLD) != 1:
        print("ERROR: original block not found exactly once (found %d) -- file differs from the 2018-01-31 revision" % text.count(OLD))
        sys.exit(1)
    with open(path + ".orig", "w", newline="") as f:
        f.write(src)
    text = text.replace(OLD, NEW)
    if "#include <set>" not in text:
        text = text.replace("#include <map>", "#include <map>\n#include <set>", 1) if "#include <map>" in text else text
    with open(path, "w", newline="") as f:
        f.write(text.replace("\n", "\r\n") if crlf else text)
    print("patched:", path, "(backup in", path + ".orig)")
