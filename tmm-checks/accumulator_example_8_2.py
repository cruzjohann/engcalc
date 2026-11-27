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


def P_acc(Cs, w, R=0.0):
    """Accumulator, optionally with a linear resistance R (s/m2) in its connection:
    flow into the vessel = h / (R + 1/(j w Cs))."""
    if R == 0.0:
        return np.array([[1, -1j * w * Cs], [0, 1]])
    return np.array([[1, -1.0 / (R + 1.0 / (1j * w * Cs))], [0, 1]])


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


def solve(f, Cs=0.0, x_acc=None, n_pts=0, R=0.0):
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
            mats.append(P_acc(Cs, w, R)); xs.append(xs[-1])
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


def peaks(Cs=0.0, x_acc=None, f=np.linspace(0.02, 2.5, 25000), R=0.0):
    hr = np.array([solve(fi, Cs, x_acc, R=R) for fi in f])
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

    # ---- Where the vessel makes it worse (0.5 m3, no losses)
    Cs = Cs_gas(0.5)
    print("\nForcing frequency: no vessel / vessel at valve / vessel at tee (0.5 m3)")
    for fr in (0.255, 0.541, 4.189 / (2 * np.pi), 1.0):
        print(f"  {fr:.3f} Hz: {solve(fr):.3f} / {solve(fr, Cs, L1 + L2):.3f} / {solve(fr, Cs, L1):.3f}")

    # ---- Higher frequencies, up to 5 Hz (0.5 m3)
    f5 = np.linspace(0.02, 5.0, 50000)
    print("\nPeaks up to 5 Hz (all reach 0.8 in this lossless model), and typical hr between peaks, 2-5 Hz")
    for name, cs_, xa_ in (("no vessel", 0.0, None), ("at the tee", Cs, L1), ("at the valve", Cs, L1 + L2)):
        hr5 = np.array([solve(fi, cs_, xa_) for fi in f5])
        pk = [f5[j] for j in range(1, len(f5) - 1) if hr5[j] > hr5[j - 1] and hr5[j] > hr5[j + 1] and hr5[j] > 0.05]
        print(f"  {name:12s}: peaks {np.round(pk, 2)}, median hr 2-5 Hz {np.median(hr5[f5 >= 2]):.3f}")

    # ---- Losses in the vessel connection: linear resistance R, 2 m3 at the valve
    Cs2 = Cs_gas(2.0)
    fg = np.linspace(0.02, 1.5, 3000)
    print("\nHighest hr from 0.02 to 1.5 Hz, 2 m3 at the valve, resistance R in the connection")
    Rs = np.array([0, 10, 20, 35, 50, 75, 100, 150, 200, 300, 400, 600, 800])
    worst = []
    for R in Rs:
        hr = np.array([solve(fi, Cs2, L1 + L2, R=R) for fi in fg])
        worst.append(hr.max())
        print(f"  R = {R:4.0f} s/m2: {hr.max():.3f} at {fg[hr.argmax()]:.3f} Hz")
    fpk = peaks(Cs2, L1 + L2, f=np.linspace(0.05, 0.4, 4000))[0][0]
    print(f"  1/(w Cs) at the lowest resonance ({fpk:.3f} Hz): {1 / (2 * np.pi * fpk * Cs2):.0f} s/m2")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4), dpi=150, gridspec_kw=dict(width_ratios=[2.2, 1]))
    a1.plot(fg, [solve(fi) for fi in fg], "k--", lw=1.2, label="No vessel")
    a1.plot(fg, [solve(fi, Cs2, L1 + L2) for fi in fg], "C1-", lw=1.3, label="2 m³ at the valve, no loss")
    a1.plot(fg, [solve(fi, Cs2, L1 + L2, R=50.0) for fi in fg], "C2-", lw=2.0, label="2 m³ at the valve, R = 50 s/m²")
    a1.set_xlim(0, 1.5); a1.set_ylim(0, 0.85); a1.grid(True, lw=0.4)
    a1.set_xlabel("Forcing frequency (Hz)"); a1.set_ylabel(r"$h_r = 2|h|/H_o$ at valve")
    a1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, fontsize=8.5, frameon=False)
    a1.set_title("Example 8-2: loss in the vessel connection", fontsize=10)
    a2.plot(Rs, worst, "C2o-", lw=1.5, ms=4)
    a2.set_xlabel("Connection resistance R (s/m²)"); a2.set_ylabel("Highest hr, 0 to 1.5 Hz")
    a2.set_ylim(0, 0.85); a2.grid(True, lw=0.4); a2.set_title("There is a best size", fontsize=10)
    fig.tight_layout(); fig.savefig("accumulator_damping.png")
