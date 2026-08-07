"""Generate three explanatory GIF animations for the COSC3500 project brief.

1. heston-paths.gif     - one simulated path of the two coupled Heston processes
2. value-surface.gif    - option value surface V(S, vol) evolving backwards from expiry
3. instability.gif      - error field exploding when dt exceeds the stability limit

Small Python reference implementations only (Euler path simulator + explicit
finite-difference solver); the real deliverables are rebuilt from the C++ engine.
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, SymLogNorm
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path

OUT = Path("/Users/zer0/Documents/cosc3500/assets")
OUT.mkdir(exist_ok=True)

# ---- palette (dataviz reference instance, light mode) ----
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
BLUE = "#2a78d6"    # categorical slot 1
ORANGE = "#eb6834"  # categorical slot 2

SEQ = LinearSegmentedColormap.from_list("seq_blue", [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
    "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b",
])
DIV = LinearSegmentedColormap.from_list("div_blue_red", [
    "#0d366b", "#1c5cab", "#2a78d6", "#86b6ef", "#cde2fb",
    "#f0efec",
    "#f6c8c7", "#f2a3a2", "#e34948", "#b23535", "#7c2222",
])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.edgecolor": BASELINE,
    "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "sans-serif", "font.size": 10,
})

# ---- model parameters (the reference configuration in the brief) ----
S0, K = 5200.0, 5250.0
R, Q = 0.045, 0.013
V0, KAPPA, THETA, XI, RHO = 0.04, 1.5, 0.04, 0.35, -0.70
T = 0.25

# =====================================================================
# Animation 1: correlated paths
# =====================================================================

def sim_path(seed, n=630):
    rng = np.random.default_rng(seed)
    dt = T / n
    S = np.empty(n + 1)
    v = np.empty(n + 1)
    S[0], v[0] = S0, V0
    z1 = rng.standard_normal(n)
    z2 = RHO * z1 + np.sqrt(1 - RHO**2) * rng.standard_normal(n)
    for k in range(n):
        vp = max(v[k], 0.0)
        S[k + 1] = S[k] * np.exp((R - Q - 0.5 * vp) * dt + np.sqrt(vp * dt) * z1[k])
        v[k + 1] = v[k] + KAPPA * (THETA - vp) * dt + XI * np.sqrt(vp * dt) * z2[k]
    return S, np.maximum(v, 0.0)


def pick_seed():
    fallback, fallback_vol = 0, 0.0
    for thresh_vol, thresh_dip in [(0.30, 0.94), (0.27, 0.95), (0.24, 0.97)]:
        for seed in range(4000):
            S, v = sim_path(seed)
            volmax = np.sqrt(v.max())
            if volmax > fallback_vol:
                fallback, fallback_vol = seed, volmax
            if volmax > thresh_vol and S.min() / S0 < thresh_dip:
                return seed
    return fallback


def anim_paths():
    seed = pick_seed()
    S, v = sim_path(seed)
    vol = np.sqrt(v) * 100.0
    n = len(S) - 1
    days = np.linspace(0, 63, n + 1)  # ~63 trading days in 3 months

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 5.0), dpi=100, sharex=True)
    fig.suptitle("Heston: two coupled random processes", fontsize=13,
                 fontweight="bold", color=INK)
    fig.text(0.5, 0.925,
             "rho = −0.70 — volatility spikes and price drops arrive together",
             ha="center", fontsize=9.5, color=INK2)

    pad1 = 0.05 * (S.max() - S.min())
    pad2 = 0.08 * (vol.max() - vol.min())
    ax1.set_ylim(S.min() - pad1, S.max() + pad1)
    ax2.set_ylim(max(0, vol.min() - pad2), vol.max() + pad2)
    ax1.set_xlim(0, 63)

    ax1.axhline(S0, color=BASELINE, lw=0.8, ls="--")
    ax2.axhline(np.sqrt(THETA) * 100, color=BASELINE, lw=0.8, ls="--")
    ax1.set_ylabel("stock price S  ($)")
    ax2.set_ylabel("volatility √v  (%)")
    ax2.set_xlabel("trading days from today")
    ax1.set_title("the stock price", fontsize=10, color=BLUE, loc="left")
    ax2.set_title("its volatility — a random process of its own",
                  fontsize=10, color=ORANGE, loc="left")
    ax1.grid(axis="x", visible=False)
    ax2.grid(axis="x", visible=False)

    l1, = ax1.plot([], [], color=BLUE, lw=1.8)
    l2, = ax2.plot([], [], color=ORANGE, lw=1.8)

    nf = 70
    frames = list(range(nf)) + [nf - 1] * 15

    def update(f):
        k = int(round(f / (nf - 1) * n))
        l1.set_data(days[: k + 1], S[: k + 1])
        l2.set_data(days[: k + 1], vol[: k + 1])
        return l1, l2

    fig.tight_layout(rect=(0, 0, 1, 0.90))
    ani = FuncAnimation(fig, update, frames=frames, blit=False)
    ani.save(OUT / "heston-paths.gif", writer=PillowWriter(fps=14))
    plt.close(fig)
    print(f"heston-paths.gif done (seed {seed}, "
          f"max vol {np.sqrt(v.max()) * 100:.0f}%, min S {S.min():.0f})")


# =====================================================================
# Explicit finite-difference solver (shared by animations 2 and 3)
# =====================================================================

NS, NV = 121, 61
S_MAX, V_MIN, V_MAX = 10400.0, 0.0, 0.5
S_AX = np.linspace(0.0, S_MAX, NS)
V_AX = np.linspace(V_MIN, V_MAX, NV)
DS = S_AX[1] - S_AX[0]
DV = V_AX[1] - V_AX[0]
DT_LIMIT = 1.0 / (2 * (0.5 * V_MAX * S_MAX**2) / DS**2 + 2 * (0.5 * XI**2 * V_MAX) / DV**2)

SI = S_AX[None, 1:-1]
VI = V_AX[1:-1, None]


def payoff():
    return np.tile(np.maximum(S_AX - K, 0.0), (NV, 1))


def step(V, dt, tau_next):
    """One explicit timestep. V has shape (NV, NS); tau_next is the new time to expiry."""
    V_S = (V[1:-1, 2:] - V[1:-1, :-2]) / (2 * DS)
    V_SS = (V[1:-1, 2:] - 2 * V[1:-1, 1:-1] + V[1:-1, :-2]) / DS**2
    V_v = (V[2:, 1:-1] - V[:-2, 1:-1]) / (2 * DV)
    V_vv = (V[2:, 1:-1] - 2 * V[1:-1, 1:-1] + V[:-2, 1:-1]) / DV**2
    V_Sv = (V[2:, 2:] - V[2:, :-2] - V[:-2, 2:] + V[:-2, :-2]) / (4 * DS * DV)
    L = ((R - Q) * SI * V_S
         + KAPPA * (THETA - VI) * V_v
         + 0.5 * VI * SI**2 * V_SS
         + 0.5 * XI**2 * VI * V_vv
         + RHO * XI * VI * SI * V_Sv
         - R * V[1:-1, 1:-1])
    Vn = V.copy()
    Vn[1:-1, 1:-1] = V[1:-1, 1:-1] + dt * L
    Vn[:, 0] = 0.0
    Vn[:, -1] = S_MAX * np.exp(-Q * tau_next) - K * np.exp(-R * tau_next)
    Vn[0, :] = Vn[1, :]
    Vn[-1, :] = Vn[-2, :]
    return Vn


# =====================================================================
# Animation 2: value surface flowing backwards from expiry
# =====================================================================

def anim_surface():
    dt = 0.9 * DT_LIMIT
    n_t = int(np.ceil(T / dt))
    dt = T / n_t

    n_keep = 76
    keep = np.unique(np.linspace(0, n_t, n_keep).astype(int))
    V = payoff()
    snaps = {0: V.copy()}
    for k in range(1, n_t + 1):
        V = step(V, dt, k * dt)
        if k in set(keep):
            snaps[k] = V.copy()
    taus = [k * dt for k in keep]
    # zoom the display window around the strike so the kink smoothing reads
    js = (S_AX >= 2600) & (S_AX <= 7800)
    jv = np.sqrt(V_AX) * 100.0 <= 55.0
    states = [snaps[k][np.ix_(jv, js)] for k in keep]

    Sg, Volg = np.meshgrid(S_AX[js] / 1000.0, np.sqrt(V_AX[jv]) * 100.0)
    zmax = max(s.max() for s in states)
    norm = Normalize(0, zmax)

    fig = plt.figure(figsize=(7.2, 5.6), dpi=100)
    ax = fig.add_subplot(projection="3d")

    nf = len(states)
    order = [0] * 8 + list(range(nf)) + [nf - 1] * 16

    def update(i):
        ax.clear()
        Z = states[i]
        ax.plot_surface(Sg, Volg, Z, facecolors=SEQ(norm(Z)), shade=False,
                        rstride=1, cstride=1, linewidth=0.1,
                        edgecolor=(1, 1, 1, 0.25), antialiased=False)
        ax.set_xlabel("stock price  ($k)", color=INK2, fontsize=9)
        ax.set_ylabel("volatility √v  (%)", color=INK2, fontsize=9)
        ax.set_zlabel("option value  ($)", color=INK2, fontsize=9)
        ax.set_zlim(0, zmax)
        ax.tick_params(colors=MUTED, labelsize=8)
        for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
            pane.set_facecolor(SURFACE)
            pane.set_edgecolor(GRID)
        ax.set_title("The whole spreadsheet, filled in backwards", fontsize=12,
                     fontweight="bold", color=INK, pad=14)
        ax.text2D(0.5, 0.945, f"{taus[i] * 365:5.0f} days before expiry",
                  transform=ax.transAxes, ha="center", fontsize=10, color=INK2)
        ax.view_init(elev=24, azim=-62 + 12 * i / (nf - 1))
        return ()

    ani = FuncAnimation(fig, update, frames=order, blit=False)
    ani.save(OUT / "value-surface.gif", writer=PillowWriter(fps=12))
    plt.close(fig)
    print(f"value-surface.gif done ({n_t} solver steps, dt {dt:.2e})")


# =====================================================================
# Animation 3: instability when dt exceeds the stable limit
# =====================================================================

def anim_instability():
    dt_bad = 4.0 * DT_LIMIT
    dt_ref = dt_bad / 8.0
    n_steps = 58

    V_u = payoff()
    V_r = payoff()
    frames = []
    for k in range(n_steps):
        tau0 = k * dt_bad
        V_u = np.clip(step(V_u, dt_bad, tau0 + dt_bad), -1e15, 1e15)
        for j in range(8):
            V_r = step(V_r, dt_ref, tau0 + (j + 1) * dt_ref)
        frames.append(((k + 1) * dt_bad, (V_u - V_r).copy()))

    norm = SymLogNorm(linthresh=1.0, linscale=0.8, vmin=-1e12, vmax=1e12, base=10)

    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=100)
    err0 = np.clip(frames[0][1], -1e12, 1e12)
    im = ax.imshow(err0, origin="lower", aspect="auto", cmap=DIV, norm=norm,
                   extent=(0, S_MAX / 1000.0, np.sqrt(V_MIN) * 100, np.sqrt(V_MAX) * 100))
    ax.grid(False)
    ax.set_xlabel("stock price  ($k)")
    ax.set_ylabel("volatility √v  (%)")
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("V − stable reference  ($)", color=INK2, fontsize=9)
    cb.ax.tick_params(colors=MUTED, labelsize=8)
    cb.outline.set_edgecolor(BASELINE)
    fig.suptitle("Take too big a timestep and the sheet explodes", fontsize=12.5,
                 fontweight="bold", color=INK)
    sub = ax.set_title("", fontsize=9.5, color=INK2)

    order = list(range(len(frames))) + [len(frames) - 1] * 14

    def update(i):
        tau, err = frames[i]
        im.set_data(np.clip(err, -1e12, 1e12))
        m = np.abs(err).max()
        sub.set_text(f"Δt = 4× the stable limit — "
                     f"step {i + 1}, max error {m:,.1e} $")
        return im, sub

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    ani = FuncAnimation(fig, update, frames=order, blit=False)
    ani.save(OUT / "instability.gif", writer=PillowWriter(fps=10))
    plt.close(fig)
    print(f"instability.gif done (final max error {np.abs(frames[-1][1]).max():.2e})")


if __name__ == "__main__":
    import sys
    print(f"stable dt limit ≈ {DT_LIMIT:.2e} years")
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "paths"):
        anim_paths()
    if which in ("all", "surface"):
        anim_surface()
    if which in ("all", "instability"):
        anim_instability()
    for f in sorted(OUT.glob("*.gif")):
        print(f"{f.name}: {f.stat().st_size / 1e6:.1f} MB")
