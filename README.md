# FCC Coral Health Report

Longitudinal analysis of coral reef colony health at Cap Cana, Dominican Republic, based on self-collected photoquadrat surveys of 5×5 m plots repeated over time. The project tracks individual coral colonies across surveys to quantify growth, partial mortality, complete mortality, and community change, and compares two reef zones: **S1P1 (Juanillo)** and **S2P2 (Farallon)**.

Colony-level measurements are derived from photogrammetry and orthoprojection segmentation of the survey imagery (see the workflow guides in the repository root), then combined and modelled with mixed-effects models to separate the effects of colony size, species, time, and site.

## Repository layout

- **`code/`** — analysis pipeline: `combine_data.py` (merges per-site survey data), `analysis.r` (statistical models and figures), and `README.md` with the full statistical results summary.
- **`data/`** — cleaned per-site colony data (`S1P1_combined_data.csv`, `S2P2_combined_data.csv`).
- **`raw data/`** — original per-site survey exports (`S1P1/`, `S2P2/`).
- **`plots/`** — generated figures (species overview, growth and mortality trajectories, size effects, site comparisons) with `figure_descriptions.md`.
- **Workflow guides** — `Guide - Photogrametry Production ES.docx`, `Guide - Orthoprojection Segmentation.docx`.

## Analysis scope

Eight species with ≥ 5 colonies are analysed (PAST, SINT, PPOR, PSTR, SSID, OFAV, DLAB, CNAT), with *Porites astreoides* (PAST) as the reference level. The modelling covers community cover and abundance, species growth trajectories, mortality severity and risk, initial-size effects, and Juanillo-vs-Farallon site contrasts. See [`code/README.md`](code/README.md) for model specifications and results.
