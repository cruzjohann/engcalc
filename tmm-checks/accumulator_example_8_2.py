"""
Gas-charged accumulator (air chamber) in Chaudhry Example 8-2.

Point matrix for an accumulator, state {q, h} (Chaudhry Problem 8-7 asks for it):
    head is the same on both sides, and the vessel takes q_acc = j*w*Cs*h, so
    {q, h}_R = [[1, -j*w*Cs], [0, 1]] {q, h}_L
    Cs = rho*g*V/(n*p): liquid volume taken in per metre of head rise
    V = gas volume at the mean absolute pressure p, n = polytropic exponent.

Uses field_matrix() and the Example 8-2 data from chaudhry_example_8_2.py.
Johann Cruz
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import chaudhry_example_8_2 as m

rho, patm = 999.0, 101325.0
g = m.g
L1, L2 = m.PIPE1["L"], m.PIPE2["L"]
P_LINE = rho * g * m.Ho + patm                 # frictionless line, valve to atmosphere


def Cs_gas(V, p=P_LINE, n=1.0):
    return rho * g * V / (n * p)


def P_acc(Cs, w):
    return np.array([[1, -1j * w * Cs], [0, 1]])


def pipe(p, L):
    return dict(p, L=L)


def segments(x_acc):
    """Main line from reservoir (x = 0) to valve (x = 750 m) as a list of ('pipe', data) / ('tee',) / ('acc',)."""
    out = []
    for name, x0, x1, data in (("1", 0.0, L1, m.PIPE1), ("2", L1, L1 + L2, m.PIPE2)):
        if x_acc is not None and x0 < x_acc < x1:
            out += [("pipe", pipe(data, x_acc - x0)), ("acc",), ("pipe", pipe(data, x1 - x_acc))]
        else:
            out.append(("pipe", data))
        if name == "1":
            out.append(("tee",))
            if x_acc == L1:
                out.append(("acc",))
    if x_acc == L1 + L2:
        out.append(("acc",))
    return out


def solve(f, Cs=0.0, x_acc=None, n_pts=0):
    """Head amplitude at the valve (hr = 2|h|/Ho); optionally |h| along the main line."""
    w = 2 * np.pi * f
    Fb = m.field_matrix(m.BRANCH, w)
    mats, xs = [], [0.0]
    for s in segments(x_acc if Cs > 0 else None):
        if s[0] == "pipe":
            mats.append(m.field_matrix(s[1], w)); xs.append(xs[-1] + s[1]["L"])
        elif s[0] == "tee":
            mats.append(np.array([[1, Fb[0, 1] / Fb[0, 0]], [0, 1]])); xs.append(xs[-1])
        else:
            mats.append(P_acc(Cs, w)); xs.append(xs[-1])
    U = np.eye(2)
    for M in mats:
        U = M @ U
    q1 = -(2 * m.Ho * m.k / m.tau0) / (U[1, 0] - 2 * m.Ho / m.Qo * U[0, 0])   # Eq. 8-133
    hr = 2 * abs(U[1, 0] * q1) / m.Ho
    if not n_pts:
        return hr
    # head amplitude along the line (for mode shapes)
    x_out, h_out, z, xpos = [], [], np.array([q1, 0]), 0.0
    for s, M in zip(segments(x_acc if Cs > 0 else None), mats):
        if s[0] == "pipe":
            for xi in np.linspace(0, s[1]["L"], n_pts)[1:]:
                zi = m.field_matrix(pipe(s[1], xi), w) @ z
                x_out.append(xpos + xi); h_out.append(abs(zi[1]))
            xpos += s[1]["L"]
        z = M @ z
    return hr, np.array(x_out), np.array(h_out)


def peaks(Cs=0.0, x_acc=None, f=np.linspace(0.02, 2.5, 25000)):
    hr = np.array([solve(fi, Cs, x_acc) for fi in f])
    i = [j for j in range(1, len(f) - 1) if hr[j] > hr[j - 1] and hr[j] > hr[j + 1] and hr[j] > 0.05]
    return f[i], hr[i]


if __name__ == "__main__":
    # ---- Check 1: Cs -> 0 gives the book case
    print("Check 1, wr = 2.0: book hand calc 0.170, no accumulator", round(solve(4.189 / (2 * np.pi)), 4),
          ", Cs = 1e-12 at valve", round(solve(4.189 / (2 * np.pi), 1e-12, 750.0), 4))

    # ---- Check 2: single pipe, reservoir -> pipe -> accumulator at a closed end
    p = m.PIPE1; A = np.pi * p["D"] ** 2 / 4; C = p["a"] / (g * A); Cp = g * A * p["L"] / p["a"] ** 2
    print(f"\nCheck 2: pipe L = {p['L']:.0f} m, D = {p['D']} m, a = {p['a']:.0f} m/s; pipe compliance gAL/a^2 = {Cp:.4f} m2")
    print(f"limit: Cs -> 0 gives the closed-end pipe, a/4L = {p['a']/4/p['L']:.3f} Hz")
    rows = []
    for V in (0.05, 0.5, 2.0, 10.0, 100.0):
        Cs = Cs_gas(V)
        ws = np.linspace(1e-3, np.pi * p["a"] / p["L"] * 0.999, 200000)
        U00 = np.array([(P_acc(Cs, w) @ m.field_matrix(p, w))[0, 0].real for w in ws])  # closed end: q = 0
        wn = ws[np.where(np.sign(U00[:-1]) != np.sign(U00[1:]))[0][0]]
        f_tmm = wn / 2 / np.pi
        f_lump = np.sqrt(g * A / (p["L"] * (Cs + Cp / 3))) / 2 / np.pi       # water column on a gas spring
        rows.append((V, Cs, f_tmm, f_lump))
        print(f"V = {V:6.2f} m3, Cs = {Cs:.4f} m2: TMM {f_tmm:.3f} Hz, lumped {f_lump:.3f} Hz, "
              f"diff {100 * (f_lump / f_tmm - 1):+.1f}%")

    # ---- Example 8-2 with an accumulator
    print(f"\nExample 8-2, line pressure {P_LINE / 1e5:.1f} bar abs, isothermal gas")
    f0, h0 = peaks()
    print("no accumulator: peaks", np.round(f0, 3), "heights", np.round(h0, 2))
    for V in (0.5, 2.0):
        Cs = Cs_gas(V)
        for where, x in (("tee", L1), ("valve", L1 + L2)):
            fp, hp = peaks(Cs, x)
            print(f"V = {V} m3 (Cs = {Cs:.4f}) at {where:5s}: peaks {np.round(fp, 3)} heights {np.round(hp, 2)}, "
                  f"hr at 0.667 Hz = {solve(4.189 / (2 * np.pi), Cs, x):.3f}, hr at 1.000 Hz = {solve(1.0, Cs, x):.3f}")

    # ---- Figure: frequency response with 0.5 m3 at the tee and at the valve
    Cs = Cs_gas(0.5)
    f = np.linspace(0.02, 2.0, 20000)
    fig, ax = plt.subplots(figsize=(8, 4.6), dpi=150)
    for x, st, lab in ((None, "k-", "No accumulator (book)"), (L1, "C0--", "0.5 m³ gas at the tee"),
                       (L1 + L2, "C3-", "0.5 m³ gas at the valve")):
        ax.plot(f, [solve(fi, Cs if x else 0.0, x) for fi in f], st, lw=1.3, label=lab)
    ax.set_xlim(0, 2.0); ax.set_ylim(0, 0.85)
    ax.set_xlabel("Forcing frequency (Hz)"); ax.set_ylabel(r"$h_r = 2|h|/H_o$ at valve")
    ax.set_title(f"Example 8-2 with an accumulator, line pressure {P_LINE / 1e5:.1f} bar abs")
    ax.grid(True, lw=0.4)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig("accumulator_response.png")

    # ---- Figure: mode shapes and resonance against accumulator position
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7.2), dpi=150, sharex=True)
    for fr, st, lab in ((0.330, "k-", "Mode 1, 0.330 Hz"), (1.0, "C0--", "Mode 2, 1.000 Hz")):
        _, xx, hh = solve(fr, n_pts=200)
        a1.plot(xx, hh / hh.max(), st, lw=1.4, label=lab)
    a1.axvline(L1, color="0.5", lw=0.8, ls=":"); a1.text(L1 + 8, 0.92, "tee", fontsize=8, color="0.35")
    a1.set_ylabel("|h| / max (no accumulator)"); a1.set_ylim(0, 1.05); a1.grid(True, lw=0.4)
    a1.set_title("Head amplitude along the main line, no accumulator", fontsize=10)
    a1.legend(loc="lower right", fontsize=9)
    xs = np.linspace(5, L1 + L2, 76)
    fr1, fr2 = [], []
    for x in xs:
        fp, _ = peaks(Cs, x, f=np.linspace(0.05, 1.3, 6000))
        fr1.append(fp[0]); fr2.append(fp[1])
    a2.plot(xs, fr1, "k-", lw=1.4, label="First resonance")
    a2.plot(xs, fr2, "C0--", lw=1.4, label="Second resonance")
    a2.axhline(0.330, color="k", lw=0.6, ls=":"); a2.axhline(1.0, color="C0", lw=0.6, ls=":")
    a2.axvline(L1, color="0.5", lw=0.8, ls=":")
    a2.set_xlabel("Accumulator position from reservoir (m)"); a2.set_ylabel("Resonant frequency (Hz)")
    a2.set_xlim(0, L1 + L2); a2.grid(True, lw=0.4)
    a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=9, frameon=False)
    a2.set_title("0.5 m³ of gas at 10.8 bar abs placed at each point (dotted lines: no accumulator)", fontsize=10)
    fig.tight_layout(); fig.savefig("accumulator_placement.png")
    np.savetxt("placement.csv", np.c_[xs, fr1, fr2], delimiter=",", fmt="%.4f", header="x,f1,f2")
