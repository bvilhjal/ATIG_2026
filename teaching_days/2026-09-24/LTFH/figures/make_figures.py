"""Regenerate the three LT-FH exercise figures into this folder.

thresholds.png     Part B: liability thresholds T40 and T70, whole population
                   and the disease-free-at-40 probands.
observed_scale.png Part E: a liability split into one 0/1 observed outcome.
h2_by_design.png   Part E answer: observed-scale h2 against sample case
                   proportion P, for liability h2 = 0.5 and K = 0.10.

Run with any Python that has NumPy and Matplotlib:
    python figures/make_figures.py
"""
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

BLUE, ORANGE, FILL, GREY, RED = "#2a78d6", "#eb6834", "#dbe9f7", "#4a4a46", "#c0392b"

# the exercise's generating model
CIP_K, CIP_MID, CIP_SLOPE = 0.10, 60.0, 1.0 / 8.0
AGE = np.arange(0, 121, 1.0)
CIP = CIP_K / (1.0 + np.exp((CIP_MID - AGE) * CIP_SLOPE))
T40 = norm.ppf(1.0 - np.interp(40.0, AGE, CIP))
T70 = norm.ppf(1.0 - np.interp(70.0, AGE, CIP))
K70 = float(np.interp(70.0, AGE, CIP))
K40 = float(np.interp(40.0, AGE, CIP))


def style(ax, xlab="liability  (higher = earlier onset)"):
    ax.set_xlabel(xlab, fontsize=9, color=GREY)
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GREY)
    ax.tick_params(colors=GREY, labelsize=8)


x = np.linspace(-3.9, 3.9, 500)
pdf = norm.pdf(x)

# ------------------------------------------------ threshold figure (Part B)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.6, 2.9))

a1.fill_between(x, pdf, color=FILL)
a1.plot(x, pdf, color=BLUE, lw=1.8)
a1.fill_between(x, pdf, where=x >= T70, color=RED, alpha=0.55, lw=0)
a1.axvline(T40, color=GREY, lw=1.0, ls="--")
a1.axvline(T70, color=GREY, lw=1.0)
a1.annotate(f"T40\n({K40:.1%} by 40)", (T40, 0.33), xytext=(T40 - 0.08, 0.475),
            ha="right", va="top", fontsize=8.5, color=GREY,
            arrowprops=dict(arrowstyle="-", color=GREY, lw=0.7))
a1.annotate(f"T70\n({K70:.0%} by 70)", (T70, 0.33), xytext=(T70 + 0.08, 0.475),
            ha="left", va="top", fontsize=8.5, color=GREY,
            arrowprops=dict(arrowstyle="-", color=GREY, lw=0.7))
a1.set_ylim(0, 0.62)
a1.text(2.75, 0.10, "cases", fontsize=8.5, color=RED, ha="center")
a1.set_title("whole population", fontsize=10)
style(a1)

# probands still undiagnosed at 40: their liability is truncated at T40
keep = x <= T40
xk, pdfk = x[keep], pdf[keep] / norm.cdf(T40)
a2.fill_between(xk, pdfk, color=FILL)
a2.plot(xk, pdfk, color=BLUE, lw=1.8)
a2.fill_between(xk, pdfk, where=xk >= T70, color=RED, alpha=0.55, lw=0)
a2.set_xlim(a1.get_xlim())
a2.axvline(T40, color=GREY, lw=1.0, ls="--")
a2.axvline(T70, color=GREY, lw=1.0)
a2.text(T70, 0.46, "T70", ha="center", fontsize=8.5, color=GREY)
a2.text(2.35, 0.07, "incident\ncases", fontsize=8.5, color=RED, ha="center")
a2.set_title("probands disease-free at 40", fontsize=10)
a2.annotate("", xy=(T40, 0.38), xytext=(T40 + 0.9, 0.38),
            arrowprops=dict(arrowstyle="->", color=GREY, lw=0.9))
a2.text(T40 + 0.95, 0.38, "everyone is\nleft of T40", fontsize=8, color=GREY, va="center")
style(a2)

fig.suptitle("One liability, two thresholds: risk is the interval between them",
             fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig("thresholds.png", dpi=200)
plt.close(fig)

# ------------------------------------------- observed-scale figure (Part E)
K, T = CIP_K, norm.ppf(1.0 - CIP_K)
fig, ax = plt.subplots(figsize=(7.4, 2.9))
ax.fill_between(x, pdf, color=FILL)
ax.plot(x, pdf, color=BLUE, lw=1.8)
ax.fill_between(x, pdf, where=x >= T, color=RED, alpha=0.55, lw=0)
ax.axvline(T, color=GREY, lw=1.0)
ax.text(T, 0.45, f"T = Φ⁻¹(1 − K)\n(K = {K:.0%})", ha="center", fontsize=8.5, color=GREY)
ax.text(-2.35, 0.26, "controls:  y = 0", fontsize=9.5, color=BLUE, ha="center")
ax.text(2.6, 0.26, "cases:  y = 1", fontsize=9.5, color=RED, ha="center")
ax.set_ylim(0, 0.5)
ax.set_title("The 0/1 observed scale: the same y for everyone on one side of T",
             fontsize=11)
style(ax, xlab="liability")
fig.tight_layout()
fig.savefig("observed_scale.png", dpi=200)
plt.close(fig)

# -------------------------------------- h2 by design curve (Part E, Q12)
z = norm.pdf(T)
P = np.linspace(0.01, 0.99, 400)
h2_obs = 0.5 * (z * z) * P * (1 - P) / (K * (1 - K)) ** 2
fig, ax = plt.subplots(figsize=(7.0, 3.1))
ax.plot(P, h2_obs, color=BLUE, lw=2)
ax.axhline(0.5, color=GREY, lw=1.0, ls="--")
ax.text(0.015, 0.505, "liability-scale h² = 0.5 (one number, any design)",
        fontsize=8.5, color=GREY, va="bottom")
for p, lbl in ((K, "population\nP = K = 0.10"), (0.2, "1:4 study\nP = 0.20"),
               (0.5, "1:1 study\nP = 0.50")):
    v = 0.5 * z * z * p * (1 - p) / (K * (1 - K)) ** 2
    ax.plot([p], [v], "o", ms=5, color=ORANGE)
    ax.annotate(f"{lbl}\nh²_obs = {v:.2f}", (p, v), textcoords="offset points",
                xytext=(6, -26 if p == 0.5 else 8), fontsize=8, color=GREY)
ax.set_xlabel("sample case proportion  P", fontsize=9, color=GREY)
ax.set_ylabel("observed-scale h²", fontsize=9, color=GREY)
ax.set_ylim(0, 0.56)
ax.set_title("The same disease reports a different observed-scale h² in every design",
             fontsize=11)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.tick_params(colors=GREY, labelsize=8)
ax.spines["bottom"].set_color(GREY)
ax.spines["left"].set_color(GREY)
fig.tight_layout()
fig.savefig("h2_by_design.png", dpi=200)
plt.close(fig)
print("wrote thresholds.png, observed_scale.png, h2_by_design.png")
