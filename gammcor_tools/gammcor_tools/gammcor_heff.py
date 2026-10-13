import sys

import numpy as np


def read(files):
    H = {}
    for name in files:
        for line in open(name):
            p = line.split()
            if len(p) == 4 and p[0] == "H_mn":
                H[(p[1], p[2])] = float(p[3].replace("D", "E"))
    return H


def order(label):
    return tuple(int(x) for x in label.split("."))


def show(title, A, states):
    print(f"\n{title}")
    print(" " * 8 + "".join(f"{s:>20}" for s in states))
    for s, row in zip(states, A):
        print(f"{s:>8}" + "".join(f"{x:20.10f}" for x in row))


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: gammcor_heff output_1 output_2 ...")

    H = read(sys.argv[1:])
    states = sorted({s for pair in H for s in pair}, key=order)
    missing = [f"H({j},{i})" for j in states for i in states if (j, i) not in H]
    if missing:
        sys.exit("missing H_mn elements: " + ", ".join(missing))


    M = np.array([[H[(j, i)] for i in states] for j in states])
    S = 0.5 * (M + M.T)
    show("H_eff (rows J, columns I)", M, states)
    show("(H_eff + H_eff^T)/2", S, states)

    print("\nEigenvalues of (H_eff + H_eff^T)/2, without ENuc")
    for k, e in enumerate(np.linalg.eigvalsh(S), 1):
        print(f"  S{k} = {e:20.10f}")


if __name__ == "__main__":
    main()
