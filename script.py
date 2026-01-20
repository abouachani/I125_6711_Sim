import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

# PARAMÈTRES SOURCE 6711
L = 0.28  # longueur active en cm
R0_CM = 1.0
THETA0 = math.pi / 2


def parse_bnn_lis(path: Path):
    header_re = {
        "r": re.compile(
            r"R coordinate: from\s+([0-9Ee+\-.]+)\s+to\s+([0-9Ee+\-.]+)\s+cm,\s+(\d+)\s+bins"
        ),
        "p": re.compile(
            r"P coordinate: from\s+([0-9Ee+\-.]+)\s+to\s+([0-9Ee+\-.]+)\s+rad,\s+(\d+)\s+bins"
        ),
        "z": re.compile(
            r"Z coordinate: from\s+([0-9Ee+\-.]+)\s+to\s+([0-9Ee+\-.]+)\s+cm,\s+(\d+)\s+bins"
        ),
    }

    r_min = r_max = p_min = p_max = z_min = z_max = None
    nr = np_ = nz = None
    data_values = []
    data_started = False

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not data_started:
                for key, regex in header_re.items():
                    match = regex.search(line)
                    if match:
                        start, end, bins = match.groups()
                        start = float(start)
                        end = float(end)
                        bins = int(bins)
                        if key == "r":
                            r_min, r_max, nr = start, end, bins
                        elif key == "p":
                            p_min, p_max, np_ = start, end, bins
                        elif key == "z":
                            z_min, z_max, nz = start, end, bins
                if "Data follow" in line:
                    data_started = True
                continue

            data_values.extend(float(val) for val in line.split())

    if None in (r_min, r_max, p_min, p_max, z_min, z_max, nr, np_, nz):
        raise ValueError(f"Missing header metadata in {path}")

    expected = nr * np_ * nz
    if len(data_values) < expected:
        raise ValueError(
            f"Not enough data values in {path}. Expected {expected}, got {len(data_values)}."
        )

    data = np.array(data_values[:expected], dtype=float)
    data = data.reshape((nr, np_, nz))

    r_edges = np.linspace(r_min, r_max, nr + 1)
    p_edges = np.linspace(p_min, p_max, np_ + 1)
    z_edges = np.linspace(z_min, z_max, nz + 1)

    r_centers = (r_edges[:-1] + r_edges[1:]) / 2
    p_centers = (p_edges[:-1] + p_edges[1:]) / 2
    z_centers = (z_edges[:-1] + z_edges[1:]) / 2

    return {
        "data": data,
        "r_centers": r_centers,
        "p_centers": p_centers,
        "z_centers": z_centers,
    }


def G_L(r, theta):
    return (1 / (L * r * np.sin(theta))) * (
        np.arctan((L / 2 + r * np.cos(theta)) / (r * np.sin(theta)))
        - np.arctan((-L / 2 + r * np.cos(theta)) / (r * np.sin(theta)))
    )


def interp_dose(r_centers, z_centers, data_rz, r, z):
    if r < r_centers[0] or r > r_centers[-1] or z < z_centers[0] or z > z_centers[-1]:
        return np.nan

    z_index = np.searchsorted(z_centers, z)
    if z_index == 0:
        z1 = z2 = 0
    elif z_index == len(z_centers):
        z1 = z2 = len(z_centers) - 1
    else:
        z1 = z_index - 1
        z2 = z_index

    dose_z1 = np.interp(r, r_centers, data_rz[:, z1])
    dose_z2 = np.interp(r, r_centers, data_rz[:, z2])

    if z1 == z2:
        return dose_z1

    return np.interp(z, [z_centers[z1], z_centers[z2]], [dose_z1, dose_z2])


def dose_at_r_theta(r_centers, z_centers, data_rz, r, theta):
    r_cyl = r * math.sin(theta)
    z = r * math.cos(theta)
    return interp_dose(r_centers, z_centers, data_rz, r_cyl, z)


def compute_tg43_parameters(sim_data, r0=R0_CM, theta0=THETA0, r_max=10.0, theta_step=5.0):
    r_centers = sim_data["r_centers"]
    z_centers = sim_data["z_centers"]
    data_rz = sim_data["data"].mean(axis=1)

    d_ref = dose_at_r_theta(r_centers, z_centers, data_rz, r0, theta0)
    if np.isnan(d_ref) or d_ref == 0:
        raise ValueError("Reference dose is invalid or zero.")

    g_ref = G_L(r0, theta0)

    r_vals = r_centers[r_centers <= r_max]
    gL_rows = []
    for r in r_vals:
        dose = dose_at_r_theta(r_centers, z_centers, data_rz, r, theta0)
        if np.isnan(dose) or dose == 0:
            continue
        gL = (dose / d_ref) * (g_ref / G_L(r, theta0))
        gL_rows.append({"r_cm": r, "gL": gL})

    gL_df = pd.DataFrame(gL_rows)

    theta_deg_vals = np.arange(theta_step, 180.0, theta_step)
    f_rows = []
    for theta_deg in theta_deg_vals:
        theta = math.radians(theta_deg)
        dose = dose_at_r_theta(r_centers, z_centers, data_rz, r0, theta)
        if np.isnan(dose) or dose == 0:
            continue
        f_val = (dose / d_ref) * (g_ref / G_L(r0, theta))
        f_rows.append({"theta_deg": theta_deg, "F": f_val})

    f_df = pd.DataFrame(f_rows)

    return {
        "dose_ref": d_ref,
        "gL": gL_df,
        "F": f_df,
    }


def compare_with_experimental(sim_df, exp_path, on_cols):
    exp_df = pd.read_csv(exp_path)
    merged = sim_df.merge(exp_df, on=on_cols, suffixes=("_sim", "_exp"))
    if merged.empty:
        raise ValueError(f"No overlapping rows to compare in {exp_path}")
    for col in [col for col in merged.columns if col.endswith("_sim")]:
        base_col = col.replace("_sim", "")
        exp_col = f"{base_col}_exp"
        if exp_col in merged:
            merged[f"{base_col}_diff_pct"] = (merged[col] / merged[exp_col] - 1) * 100
    return merged


def write_outputs(output_dir, label, results, exp_gL=None, exp_F=None):
    output_dir.mkdir(parents=True, exist_ok=True)

    gL_path = output_dir / f"tg43_gL_{label}.csv"
    F_path = output_dir / f"tg43_F_{label}.csv"

    results["gL"].to_csv(gL_path, index=False)
    results["F"].to_csv(F_path, index=False)

    comparisons = {}
    if exp_gL:
        comparisons["gL"] = compare_with_experimental(results["gL"], exp_gL, ["r_cm"])
        comparisons["gL"].to_csv(output_dir / f"tg43_gL_{label}_comparison.csv", index=False)
    if exp_F:
        comparisons["F"] = compare_with_experimental(results["F"], exp_F, ["theta_deg"])
        comparisons["F"].to_csv(output_dir / f"tg43_F_{label}_comparison.csv", index=False)

    return gL_path, F_path, comparisons


def main():
    parser = argparse.ArgumentParser(
        description="Calcule les paramètres TG-43 à partir des fichiers .bnn.lis."
    )
    parser.add_argument(
        "--water",
        type=Path,
        default=Path("I125_6711_51_eau_sim.bnn.lis"),
        help="Fichier .bnn.lis dans le fantôme d'eau.",
    )
    parser.add_argument(
        "--air",
        type=Path,
        default=Path("I125_6711_51_air_sim.bnn.lis"),
        help="Fichier .bnn.lis dans le fantôme d'air.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("tg43_outputs"),
        help="Répertoire de sortie pour les CSV.",
    )
    parser.add_argument(
        "--r-max",
        type=float,
        default=10.0,
        help="Rayon maximal (cm) pour gL(r).",
    )
    parser.add_argument(
        "--theta-step",
        type=float,
        default=5.0,
        help="Pas angulaire (deg) pour F(r,theta) à r=1 cm.",
    )
    parser.add_argument(
        "--exp-gL",
        type=Path,
        help="CSV expérimental avec colonnes r_cm et gL pour comparaison.",
    )
    parser.add_argument(
        "--exp-F",
        type=Path,
        help="CSV expérimental avec colonnes theta_deg et F pour comparaison.",
    )
    args = parser.parse_args()

    water_data = parse_bnn_lis(args.water)
    air_data = parse_bnn_lis(args.air)

    water_results = compute_tg43_parameters(
        water_data, r_max=args.r_max, theta_step=args.theta_step
    )
    air_results = compute_tg43_parameters(
        air_data, r_max=args.r_max, theta_step=args.theta_step
    )

    sk_proxy = air_results["dose_ref"] * (R0_CM**2)
    lambda_proxy = water_results["dose_ref"] / sk_proxy if sk_proxy != 0 else np.nan

    gL_path, f_path, comparisons = write_outputs(
        args.output_dir, "water", water_results, exp_gL=args.exp_gL, exp_F=args.exp_F
    )
    write_outputs(args.output_dir, "air", air_results)

    summary = [
        f"Dose ref eau (r=1 cm, theta=90°): {water_results['dose_ref']:.6e}",
        f"Dose ref air (r=1 cm, theta=90°): {air_results['dose_ref']:.6e}",
        f"S_k proxy (air @1cm * r^2): {sk_proxy:.6e}",
        f"Lambda proxy: {lambda_proxy:.6e}",
        f"gL CSV eau: {gL_path}",
        f"F CSV eau: {f_path}",
    ]

    if comparisons:
        summary.append("Comparaisons expérimentales générées:")
        if "gL" in comparisons:
            summary.append("- gL: tg43_gL_water_comparison.csv")
        if "F" in comparisons:
            summary.append("- F: tg43_F_water_comparison.csv")

    summary_path = args.output_dir / "tg43_summary.txt"
    summary_path.write_text("\n".join(summary) + "\n", encoding="utf-8")

    print("\n".join(summary))


if __name__ == "__main__":
    main()
