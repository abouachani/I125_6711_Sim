import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# PARAMÈTRES SOURCE 6711
L = 0.28  # longueur active en cm
r0 = 1.0
theta0 = np.pi / 2

# LECTURE DU FICHIER FLUKA
file = "I125_6711_51_plot.dat"

data = pd.read_csv(
    file,
    delim_whitespace=True,
    names=["theta", "r", "dose", "err"]
)

# FONCTION GÉOMÉTRIQUE TG-43
def G_L(r, theta):
    return (1 / (L * r * np.sin(theta))) * (
        np.arctan((L/2 + r*np.cos(theta)) / (r*np.sin(theta))) -
        np.arctan((-L/2 + r*np.cos(theta)) / (r*np.sin(theta)))
    )

# POINT DE RÉFÉRENCE
# Trouver la ligne correspondant à (r0, theta0)
ref_df = data[
    (np.isclose(data["r"], r0, atol=1e-3)) &
    (np.isclose(np.abs(data["theta"]), theta0, atol=1e-3))
]

if ref_df.empty:
    # Fallback: find rows at the reference radius and pick the theta closest to theta0
    ref_by_r = data[np.isclose(data["r"], r0, atol=1e-3)]
    if ref_by_r.empty:
        raise ValueError(
            f"Reference radius not found: r={r0}.\n"
            f"Check input file '{file}' and adjust tolerances or verify data."
        )
    # pick the row whose |abs(theta) - theta0| is minimal
    closest_idx = (np.abs(np.abs(ref_by_r["theta"]) - theta0)).idxmin()
    ref = ref_by_r.loc[closest_idx]
    print(f"Warning: exact theta={theta0} not found at r={r0}; using closest theta={ref['theta']}")
else:
    ref = ref_df.iloc[0]

D_ref = ref["dose"]
G_ref = G_L(r0, theta0)

# CALCUL gL(r)
gL = []

for r in sorted(data["r"].unique()):
    row = data[
        (np.isclose(data["r"], r, atol=1e-3)) &
        (np.isclose(np.abs(data["theta"]), theta0, atol=1e-3))
    ]

    if len(row) == 0:
        continue

    D = row["dose"].values[0]
    g = (D / D_ref) * (G_ref / G_L(r, theta0))
    gL.append([r, g])

gL = pd.DataFrame(gL, columns=["r_cm", "gL"])

# CALCUL F(r,θ) POUR r FIXÉ
r_fixed = 1.0
subset = data[np.isclose(data["r"], r_fixed, atol=1e-3)]

F = []
for _, row in subset.iterrows():
    theta = row["theta"]
    D = row["dose"]
    F_theta = (D / D_ref) * (G_ref / G_L(r_fixed, theta))
    F.append([theta * 180 / np.pi, F_theta])

F = pd.DataFrame(F, columns=["theta_deg", "F"])

# TRACÉS
plt.figure()
plt.plot(gL["r_cm"], gL["gL"], "o-", label="FLUKA")
plt.xlabel("r (cm)")
plt.ylabel("gL(r)")
plt.title("Fonction radiale gL(r) – I-125 6711")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.plot(F["theta_deg"], F["F"], "o-")
plt.xlabel("θ (deg)")
plt.ylabel("F(r,θ)")
plt.title("Fonction anisotropique F(r,θ) à r = 1 cm")
plt.grid()
plt.show()
