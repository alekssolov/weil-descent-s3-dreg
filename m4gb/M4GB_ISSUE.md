# Gebauer–Möller F-criterion in `update()`: new critical pairs with equal lcm eliminate each other

**Component:** `src/m4gb.hpp`, `m4gb<...>::update()`, the loop that filters `newCP_good`
(lines 1634–1644 in the 2018-01-31 revision).

## Summary

When a new basis polynomial `p` is inserted, `update()` builds all new critical pairs
`(p, g)` and then applies a divisibility criterion:

```cpp
for (auto it = newCP_good.begin(); it != newCP_good.end(); ++it)
{
    if (!it->second) continue;
    for (auto it2 = newCP_good.begin(); it2 != newCP_good.end() && it2->first.lcmdeg <= it->first.lcmdeg; ++it2)
        if (it != it2 && it2->first.intlcm | it->first.intlcm)   // '|' is non-strict divisibility
        {
            it->second = false;
            break;
        }
}
```

Because `operator|` is non-strict, two *good* new pairs `A = (p, g_i)` and `B = (p, g_j)` with
`lcm(A) == lcm(B)` eliminate each other: `A` is marked bad because `lcm(B) | lcm(A)`, and `B` is
marked bad because `lcm(A) | lcm(B)` (the goodness of `it2` is not consulted). Both pairs are
discarded.

The Gebauer–Möller criteria say something different:

* (M) drop `(p, g_i)` if some other new pair has an lcm that **properly** divides `lcm(p, g_i)`;
* (F) among new pairs with **equal** lcm keep **exactly one** — none only if one of them has
  coprime leading monomials (product criterion).

Discarding *all* pairs of an lcm class removes a generator of the leading-term syzygy module.
In practice M4GB still terminates with a correct Gröbner basis on the systems we tested — the
lost S-polynomials are recovered later through pairs of higher lcm degree — but the observable
effects are (a) the maximum lcm degree reached is inflated by one on systems where the missing
S-polynomials carry the degree falls, and (b) extra work at the higher degree. Strictly
speaking the Gebauer–Möller termination argument no longer applies.

## Minimal reproducer

`wd_n16.in` (attached): 16 quadratic polynomials + 16 field equations `Xi*Xi+Xi` over GF(2)
in 16 variables (Weil descent of Semaev's third summation polynomial with linear constraints,
the setting of Huang–Kosters–Yeo, CRYPTO 2015, Sec. 5.2). The system has 4 solutions; its
reduced Gröbner basis consists of 14 linear polynomials.

```
make bin/solver_m4gb.exe MAXVARS=16 FIELDSIZE=2
./bin/solver_m4gb -s -i wd_n16.in 2>&1 | grep -E "lcmdegree|#n:"
```

Before the fix:

```
P lcmdegree=3   P #n:76
P lcmdegree=3   P #n:103
P lcmdegree=4   P #n:13      <- the 13 remaining linear polynomials appear only here
P lcmdegree=3   P #n:0
```

After the fix (patch below):

```
P lcmdegree=3   P #n:76
P lcmdegree=3   P #n:13      <- found at lcm degree 3
P lcmdegree=3   P #n:0
P lcmdegree=4   P #n:0
```

The returned Gröbner basis is identical in both runs (X14, X13, X12, X11, X9+1, X8+X10+1,
X7+X15, X6, X5, X4, X3, X2+X10+1, X1+1, X0+X10, X15^2+X15, X10^2+X10).

Independent checks that lcm degree 3 is sufficient for this input: (i) the reduced Gröbner
basis lies in the degree-3 closure V_3 of the generators (in the sense of Huang–Kosters–Yeo,
Def. 1), verified by comparing with the basis of the vanishing ideal of the 4 solution points
computed with Buchberger–Möller; (ii) a plain Buchberger implementation with the normal
selection strategy (all pairs of minimal lcm degree, full interreduction, no criteria) finishes
without exceeding lcm degree 3; (iii) Magma's F4 reports step degree 3 for systems of this
family at n=16 (HKY, Table in Sec. 5.2).

Statistics: on 41 such systems (the file above, 20 random permutations/affine changes of its
variables, 10 further random instances, 10 instances without solutions) the unpatched solver
always reached lcm degree 4; the patched solver finishes at 3 in all 41 cases. On systems whose
intrinsic degree is 4 (the same family at n=18, 20) both versions report 4, so the fix does not
spuriously lower the degree.

## Proposed fix

```diff
--- src/m4gb.hpp.orig
+++ src/m4gb.hpp
@@ -1631,16 +1631,35 @@
 				else
 					newCP_good.emplace(std::move(newCP), true); // good
 			}
-			for (auto it = newCP_good.begin(); it != newCP_good.end(); ++it)
+			// Gebauer-Moeller, corrected:
+			//  (M) drop a pair if another new pair has an lcm that PROPERLY divides its lcm;
+			//  (F) among new pairs with EQUAL lcm keep exactly one -- none if one of them is coprime.
+			// (original code used non-strict divisibility, so equal-lcm pairs discarded each other)
 			{
-				if (!it->second)
-					continue;
-				for (auto it2 = newCP_good.begin(); it2 != newCP_good.end() && it2->first.lcmdeg <= it->first.lcmdeg; ++it2)
-					if (it != it2 && it2->first.intlcm | it->first.intlcm)
+				std::set<int_monomial_t> lcm_coprime, lcm_kept;
+				for (auto it = newCP_good.begin(); it != newCP_good.end(); ++it)
+					if (!it->second)
+						lcm_coprime.insert(it->first.intlcm);
+				for (auto it = newCP_good.begin(); it != newCP_good.end(); ++it)
+				{
+					if (!it->second)
+						continue;
+					if (lcm_coprime.count(it->first.intlcm) != 0)
 					{
 						it->second = false;
-						break;
+						continue;
 					}
+					for (auto it2 = newCP_good.begin(); it2 != newCP_good.end() && it2->first.lcmdeg <= it->first.lcmdeg; ++it2)
+						if (it != it2 && it2->first.intlcm != it->first.intlcm && it2->first.intlcm | it->first.intlcm)
+						{
+							it->second = false;
+							break;
+						}
+					if (!it->second)
+						continue;
+					if (!lcm_kept.insert(it->first.intlcm).second)
+						it->second = false; // an earlier pair with the same lcm is already kept
+				}
 			}
```

(`#include <set>` is needed if not already present.) `newCP_good` is ordered by `(intlcm, p2)`,
so pairs with equal lcm are contiguous and the first surviving one is kept. Builds cleanly with
`g++ -O3 -Wall -Wfatal-errors` (GCC, Cygwin).

## Additional observations

* With the fix, the solver was run on Weil-descent systems with 26–40 variables (12 runs,
  up to 6 h / 2.2 GiB each); the reported degrees agree with an independent Macaulay-matrix
  computation wherever the latter is feasible. One run (`wd_n26_s26000.in`, 26 variables)
  ended with a segmentation fault *after* the last productive step (all new polynomials had
  already been reported); the unpatched build completes the same input. The crash is
  reproducible with the patched build but Cygwin gdb could not catch it (`USETHREADS` is
  defined in `m4gb.hpp` and marked "NOT FULLY IMPLEMENTED"). We have not diagnosed it.

Attachments: `wd_n16.in`, full solver output before/after (`m4gb_full_n16_before_patch.txt`,
`m4gb_traces_session3.txt`), `m4gb_gebauer_moeller.patch`.
