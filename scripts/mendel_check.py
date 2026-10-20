#!/usr/bin/env python3
import argparse
import subprocess
from itertools import combinations

ALLELES = {0: (0,), 1: (0, 1), 2: (1,)}

def dosage(gt):
    p = gt.replace("|", "/").split("/")
    if len(p) != 2 or "." in p:
        return None
    return int(p[0]) + int(p[1])

def pair_error(a, b):
    return (a == 0 and b == 2) or (a == 2 and b == 0)

def trio_error(f, m, c):
    return c not in {x + y for x in ALLELES[f] for y in ALLELES[m]}

def count(rows, n):
    pairs = list(combinations(range(n), 2))
    trios = [(a, b, c) for c in range(n)
             for a, b in combinations([i for i in range(n) if i != c], 2)]
    sp = {k: [0, 0] for k in pairs}
    st = {k: [0, 0] for k in trios}
    for g in rows:
        for a, b in pairs:
            if g[a] is None or g[b] is None:
                continue
            sp[(a, b)][0] += 1
            sp[(a, b)][1] += pair_error(g[a], g[b])
        for a, b, c in trios:
            if g[a] is None or g[b] is None or g[c] is None:
                continue
            st[(a, b, c)][0] += 1
            st[(a, b, c)][1] += trio_error(g[a], g[b], g[c])
    return sp, st

def report(sp, st, names):
    print("## Paare (sortiert nach Fehlerrate)")
    print("probe1\tprobe2\tsites\tmendel_errors\terror_rate")
    for (a, b), (n, e) in sorted(sp.items(), key=lambda kv: kv[1][1] / max(kv[1][0], 1)):
        print(f"{names[a]}\t{names[b]}\t{n}\t{e}\t{e / n if n else float('nan'):.5f}")
    print("\n## Trios (sortiert nach Fehlerrate)")
    print("elter1\telter2\tkind\tsites\tmendel_errors\terror_rate")
    for (a, b, c), (n, e) in sorted(st.items(), key=lambda kv: kv[1][1] / max(kv[1][0], 1)):
        print(f"{names[a]}\t{names[b]}\t{names[c]}\t{n}\t{e}\t{e / n if n else float('nan'):.5f}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    a = ap.parse_args()
    names = subprocess.run(["bcftools", "query", "-l", a.vcf],
                           capture_output=True, text=True, check=True).stdout.split()
    proc = subprocess.Popen(["bcftools", "query", "-f", "[%GT\t]\n", a.vcf],
                            stdout=subprocess.PIPE, text=True)
    rows = ([dosage(g) for g in l.rstrip("\n").rstrip("\t").split("\t")] for l in proc.stdout)
    sp, st = count(rows, len(names))
    report(sp, st, names)

if __name__ == "__main__":
    main()
