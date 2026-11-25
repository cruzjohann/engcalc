"""
Frequency response of a branching pipe system by the transfer matrix method.

Reference: M. H. Chaudhry, Applied Hydraulic Transients, 3rd ed., Springer,
2014, Chapter 8, Example 8-2 and Fig. 8-17. Python mirror of
chaudhry_example_8_2.m (same equations, same results).

Johann Cruz
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

g = 9.80665                 # m/s2, standard gravity
Ho, Qo = 100.0, 0.314       # m, m3/s
k, tau0 = 0.2, 1.0          # gate oscillation amplitude, mean opening
Tth = 3.0                   # s, theoretical period 4(L1+L2)/a
wth = 2 * np.pi / Tth

PIPE1 = dict(L=500.0, a=1000.0, D=1.81)
PIPE2 = dict(L=250.0, a=1000.0, D=1.05)
BRANCH = dict(L=200.0, a=800.0, D=1.36)


def field_matrix(p, w):
    """Frictionless field matrix of a uniform pipe (Eq. 8-37)."""
    A = np.pi * p["D"] ** 2 / 4
    C = p["a"] / (g * A)
    b = p["L"] / p["a"]
    return np.array([[np.cos(b * w), -1j * np.sin(b * w) / C],
                     [-1j * C * np.sin(b * w), np.cos(b * w)]])


def response(wr):
    w = wr * wth
    F1, F2, Fb = field_matrix(PIPE1, w), field_matrix(PIPE2, w), field_matrix(BRANCH, w)
    Pbde = np.array([[1, Fb[0, 1] / Fb[0, 0]], [0, 1]])   # Eq. 8-87
    U = F2 @ Pbde @ F1                                      # Fig. 8-17b
    q1 = -(2 * Ho * k / tau0) / (U[1, 0] - 2 * Ho / Qo * U[0, 0])   # Eq. 8-133
    q3, h3 = U[0, 0] * q1, U[1, 0] * q1                     # Eqs. 8-134, 8-135
    return dict(hr=2 * abs(h3) / Ho, qr=2 * abs(q3) / Qo,
                phase_h=np.degrees(np.angle(h3)), phase_q=np.degrees(np.angle(q3)))


if __name__ == "__main__":
    r = response(2.0)
    print("Example 8-2 at wr = 2.0")
    print(f"{'Quantity':34s}{'Book':>10s}{'This code':>12s}")
    print(f"{'hr = 2|h3|/Ho':34s}{0.170:10.3f}{r['hr']:12.4f}")
    print(f"{'qr = 2|q3|/Qo':34s}{0.390:10.3f}{r['qr']:12.4f}")
    print(f"{'Phase, head vs gate (deg)':34s}{-102.29:10.2f}{r['phase_h']:12.4f}")
    print(f"{'Phase, flow vs gate (deg)':34s}{-12.29:10.2f}{r['phase_q']:12.4f}")

    wr = np.linspace(0.5, 7.0, 3001)
    res = [response(x) for x in wr]
    hr = np.array([x["hr"] for x in res]); qr = np.array([x["qr"] for x in res])

    fig, ax = plt.subplots(figsize=(8, 4.6), dpi=150)
    ax.plot(wr, hr, "k-", lw=1.4, label=r"$h_r = 2|h_3|/H_o$")
    ax.plot(wr, qr, "k--", lw=1.2, label=r"$q_r = 2|q_3|/Q_o$")
    ax.plot([2.0, 2.0], [r["hr"], r["qr"]], "o", mfc="w", mec="k", label=r"Hand-calc point, $\omega_r$ = 2.0")
    ax.set_xlim(0.5, 7); ax.set_ylim(0, 0.85)
    ax.set_xlabel(r"$\omega_r = \omega/\omega_{th}$"); ax.set_ylabel(r"$h_r,\ q_r$")
    ax.set_title("Chaudhry Example 8-2: frequency response at the oscillating valve")
    ax.grid(True, lw=0.4)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig("example_8_2_frequency_response.png")
