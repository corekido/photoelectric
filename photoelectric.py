import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2

# data
x = np.array([5.19, 5.49, 6.88, 7.41])  # nu in units of 1e14 Hz (no error)
y = np.array([0.42, 0.55, 1.04, 1.3])  # stopping voltage in V

V_error = 0.03  # V (dominated by locating the zero-current point on the analog meter)
nu_error = 0.0  # exact, so a y-only weighted fit is appropriate
sigma_y = np.full_like(y, V_error)

# weighted linear fit (weights = 1/sigma)
# cov="unscaled" -> the covariance is computed from the given sigma_y
# (not rescaled by the scatter of the points)
(slope, intercept), cov = np.polyfit(x, y, 1, w=1 / sigma_y, cov="unscaled")
slope_err = np.sqrt(cov[0, 0])
intercept_err = np.sqrt(cov[1, 1])
cov_ab = cov[0, 1]

fit_line = slope * x + intercept

# goodness of fit
dof = len(x) - 2
chi2_val = np.sum(((y - fit_line) / sigma_y) ** 2)
chi2_red = chi2_val / dof
p_value = chi2.sf(chi2_val, dof)

# 1-sigma band of the fitted line (includes slope-intercept covariance)
x_fit = np.linspace(x.min() - 0.2, x.max() + 0.2, 200)
y_fit = slope * x_fit + intercept
y_fit_err = np.sqrt(x_fit**2 * cov[0, 0] + 2 * x_fit * cov_ab + cov[1, 1])

# physical quantities: V = (h/e) * nu - phi/e
e = 1.602176634e-19  # C
h = slope * 1e-14 * e  # J s  (x is in 1e14 Hz)
h_err = slope_err * 1e-14 * e
phi = -intercept  # eV  (work function)
phi_err = intercept_err

print(f"a = {slope:.4f} +/- {slope_err:.4f}  (V per 1e14 Hz)")
print(f"b = {intercept:.3f} +/- {intercept_err:.3f}  V")
print(
    f"chi2 = {chi2_val:.2f}, dof = {dof}, chi2_red = {chi2_red:.2f}, p = {p_value:.3f}"
)
print(f"h   = ({h * 1e34:.2f} +/- {h_err * 1e34:.2f}) e-34 J s")
print(f"phi = {phi:.2f} +/- {phi_err:.2f} eV")

# style
plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
    }
)

fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

# data points with error bars
ax.errorbar(
    x,
    y,
    yerr=sigma_y,
    fmt="o",
    color="#1f77b4",
    markersize=5.5,
    markeredgecolor="black",
    markeredgewidth=0.6,
    ecolor="black",
    elinewidth=0.9,
    capsize=3,
    label="Data Points",
    zorder=3,
)
ax.plot(x, y, color="gray")

# fitted line and 1-sigma band
ax.plot(
    x_fit, y_fit, color="#d62728", linewidth=1.5, label="Weighted Linear Fit", zorder=2
)
ax.fill_between(
    x_fit,
    y_fit - y_fit_err,
    y_fit + y_fit_err,
    color="#d62728",
    alpha=0.2,
    label=r"$1\sigma$ band",
    zorder=1.5,
)

# label
ax.set_xlabel(r"$\nu \ (10^{14}\ \mathrm{Hz})$")
ax.set_ylabel(r"$V \  (\mathrm{V})$")

# page geometry
ax.tick_params(direction="in", top=True, right=True, which="both")
ax.grid(True, linestyle="--", alpha=0.4, zorder=1)

# info box
fit_info = (
    r"$y = a \cdot x + b$ " + "\n"
    rf"$a = {slope:.4f} \pm {slope_err:.4f}$" + "\n"
    rf"$b = {intercept:.3f} \pm {intercept_err:.3f}$"
)

ax.text(
    0.05,
    0.95,
    fit_info,
    transform=ax.transAxes,
    verticalalignment="top",
    bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="gray", alpha=0.9),
)

ax.legend(loc="lower right")

plt.tight_layout()

plt.savefig("photoelectric.pdf")
