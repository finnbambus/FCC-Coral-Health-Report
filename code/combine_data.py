import os
import csv
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(PROJECT_ROOT, "raw data")
DATA_DIR = os.path.join(PROJECT_ROOT, "data")

ZONE_MAP = {
    "S1P1": "Juanillo",
    "S2P2": "Farallon",
}

def make_code(species_name):
    parts = species_name.strip().split()
    if len(parts) < 2:
        return species_name[:4].upper()
    return (parts[0][0] + parts[1][:3]).upper()

# Optional: pass a site name as argument to process only that site
target_site = sys.argv[1] if len(sys.argv) > 1 else None

rows = []

for site in sorted(os.listdir(BASE)):
    site_path = os.path.join(BASE, site)
    if not os.path.isdir(site_path) or site.startswith("."):
        continue
    if target_site and site != target_site:
        continue
    zone = ZONE_MAP.get(site, site)
    for plot in sorted(os.listdir(site_path)):
        plot_path = os.path.join(site_path, plot)
        if not os.path.isdir(plot_path) or plot.startswith("."):
            continue
        for fname in sorted(os.listdir(plot_path)):
            if not fname.endswith(".csv"):
                continue
            fpath = os.path.join(plot_path, fname)
            with open(fpath, newline="", encoding="utf-8-sig") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    row = {k.strip(): v for k, v in row.items()}
                    species = row.get("TagLab Class name", "").strip()
                    rows.append([
                        zone,
                        site,
                        plot,
                        row.get("TagLab Date", "").strip(),
                        row.get("TagLab Genet Id", "").strip(),
                        species,
                        make_code(species),
                        row.get("TagLab Area", "").strip(),
                        "",  # partial_mortality
                        "",  # notes
                    ])

site_label = target_site if target_site else "all"
out_path = os.path.join(DATA_DIR, f"{site_label}_combined_data.csv")
os.makedirs(DATA_DIR, exist_ok=True)
with open(out_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["zone", "site", "plot", "date", "genet_id", "class", "code", "area", "partial_mortality", "notes"])
    writer.writerows(rows)

print(f"Written {len(rows)} rows to {out_path}")
