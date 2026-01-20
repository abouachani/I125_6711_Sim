# I125_6711_Sim

## Calcul des paramètres TG-43

Le script `script.py` extrait les paramètres TG-43 (fonction radiale gL(r) et fonction anisotropique F(r,θ)) depuis les fichiers FLUKA `.bnn.lis` pour les fantômes d’air et d’eau.

### Exécution

```bash
python script.py \
  --water I125_6711_51_eau_sim.bnn.lis \
  --air I125_6711_51_air_sim.bnn.lis \
  --output-dir tg43_outputs
```

Les sorties principales sont :

- `tg43_outputs/tg43_gL_water.csv`
- `tg43_outputs/tg43_F_water.csv`
- `tg43_outputs/tg43_gL_air.csv`
- `tg43_outputs/tg43_F_air.csv`
- `tg43_outputs/tg43_summary.txt`

### Comparaison aux données expérimentales

Vous pouvez comparer aux valeurs expérimentales du site CLRP en fournissant des CSV locaux contenant les colonnes suivantes :

- `r_cm`, `gL` (pour la comparaison de gL)
- `theta_deg`, `F` (pour la comparaison de F)

Exemple :

```bash
python script.py \
  --water I125_6711_51_eau_sim.bnn.lis \
  --air I125_6711_51_air_sim.bnn.lis \
  --output-dir tg43_outputs \
  --exp-gL OncoSeed_6711_gL.csv \
  --exp-F OncoSeed_6711_F.csv
```

Le script génère alors :

- `tg43_outputs/tg43_gL_water_comparison.csv`
- `tg43_outputs/tg43_F_water_comparison.csv`
