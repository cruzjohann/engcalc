"""
Free gas and the frequency response of Chaudhry Example 8-2.

Wave speed in a gas-liquid mixture from Chaudhry, Applied Hydraulic
Transients, 3rd ed., Section 9-5, Eq. 9-18 (small, evenly spread bubbles,
isothermal). The liquid and pipe-wall terms are taken from the gas-free
wave speed ao:  1/(rho_m am^2) = 1/(rho ao^2) + alpha/p.

Uses response() from chaudhry_example_8_2.py in this folder.
Johann Cruz
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import chaudhry_example_8_2 as m

rho, K, g, patm = 999.0, 2.19e9, 9.80665, 101325.0
R_air, T = 287.0, 293.15


def a_mix(a0, alpha, p):
    """Wave speed with free gas (Eq. 9-18). alpha = void fraction at absolute pressure p (Pa)."""
    rho_g = p / (R_air * T)
    rho_m = rho * (1 - alpha) + alpha * rho_g
    return 1 / np.sqrt(rho_m * (alpha / p + 1 / (rho * a0 ** 2)))


def set_wave_speeds(alpha, p):
    base = {"PIPE1": 1000.0, "PIPE2": 1000.0, "BRANCH": 800.0}
    for name, a0 in base.items():
        getattr(m, name)["a"] = a_mix(a0, alpha, p)
    return {name: getattr(m, name)["a"] for name in base}


if __name__ == "__main__":
    # Book statement, Section 2-6: 1 part air in 10,000 parts water halves the wave speed
    a0 = np.sqrt(K / rho)
    print(f"Rigid pipe, gas-free: {a0:.0f} m/s; 0.01% gas at 1 atm: {a_mix(a0, 1e-4, patm):.0f} m/s "
          f"(ratio {a_mix(a0, 1e-4, patm) / a0:.3f})")

    print("\nWave speed (m/s), gas-free 1000 m/s")
    for p_bar in (1.013, 3, 10.8, 30):
        print(f"{p_bar:6.3f} bar abs:", [round(a_mix(1000.0, x, p_bar * 1e5)) for x in (1e-5, 1e-4, 1e-3, 1e-2)])

    p_line = rho * g * m.Ho + patm          # frictionless line, pipe at valve level
    f = np.linspace(0.02, 1.4, 40000)
    w_hand = 4.189                            # rad/s, forcing of the book's hand calculation (0.667 Hz)
    fig, ax = plt.subplots(figsize=(8, 4.6), dpi=150)
    print(f"\nExample 8-2, line pressure {p_line / 1e5:.1f} bar abs")
    for alpha, st in [(0, "k-"), (1e-4, "C0-"), (1e-3, "C1--"), (1e-2, None)]:
        a = set_wave_speeds(alpha, p_line)
        hr = np.array([m.response(fi * 2 * np.pi / m.wth)["hr"] for fi in f])
        pk = [i for i in range(1, len(f) - 1) if hr[i] > hr[i - 1] and hr[i] > hr[i + 1] and hr[i] > 0.05]
        print(f"gas {alpha * 100:g}%: a = {[round(v) for v in a.values()]}, first resonance {f[pk[0]]:.3f} Hz, "
              f"hr at 0.667 Hz = {m.response(w_hand / m.wth)['hr']:.3f}")
        if st:
            ax.plot(f, hr, st, lw=1.3, label=f"{alpha * 100:g}% gas" if alpha else "No gas (book)")
    set_wave_speeds(0, p_line)
    ax.axvline(w_hand / (2 * np.pi), color="0.5", lw=0.8, ls=":")
    ax.set_xlim(0, 1.4); ax.set_ylim(0, 0.85)
    ax.set_xlabel("Forcing frequency (Hz)"); ax.set_ylabel(r"$h_r = 2|h|/H_o$ at valve")
    ax.set_title(f"Example 8-2 with free gas, line pressure {p_line / 1e5:.1f} bar abs")
    ax.grid(True, lw=0.4)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig("example_8_2_with_free_gas.png")
