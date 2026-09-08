"""Compose computed vector graphics and LaTeX; export Figure 28a as SVG/PNG."""
from pathlib import Path
import re
import subprocess
import tempfile
import argparse
import shutil
import os
from lxml import etree as ET

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=Path('outputs/2026-09-08/ld_score_latex'))
parser.add_argument('--offline', action='store_true', help='Use Tectonic cached packages only')
parser.add_argument('--equation-number', default=None, help='Optional current lecture equation number')
args = parser.parse_args()
OUT = args.output.resolve()
OUT.mkdir(parents=True, exist_ok=True)

def executable(env, name):
    path = os.environ.get(env) or shutil.which(name)
    if not path:
        raise SystemExit(f"Missing {name}; install it or set {env} to its executable path")
    return path

R = executable('RSCRIPT', 'Rscript')
TECTONIC = executable('TECTONIC', 'tectonic')
PDFTOCAIRO = executable('PDFTOCAIRO', 'pdftocairo')
NODE = executable('NODE', 'node')
SHARP = os.environ.get('SHARP_MODULE') or str(HERE / 'node_modules/sharp')
NAME = 'LD_score_heatmap'

FORMULAS = [
    r"\ell_i=\sum_{m=1}^{M}r_{im}^{\,2}=1+\sum_{m\ne i}r_{im}^{\,2}",
    r"\ell_4=\underbrace{\color{self}1}_{\mathrm{self}}"
    r"+\underbrace{\color{ld}2.35}_{\mathrm{other\ SNPs}}=\mathbf{3.35}",
]

with tempfile.TemporaryDirectory(prefix="atig_ld_score_") as tmp:
    tmp = Path(tmp)
    subprocess.run([R, str(HERE / "draw_ld_score.R"), str(tmp)], check=True)
    subprocess.run([PDFTOCAIRO, "-svg", str(tmp / "base.pdf"), str(tmp / "base.svg")], check=True)
    for name in ("population_genotype_correlation.csv", "population_squared_correlation.csv",
                 "ld_scores_and_contributions.csv", "model.json"):
        (OUT / name).write_bytes((tmp / name).read_bytes())
    tex = (r"\documentclass{article}\usepackage{amsmath,amssymb,xcolor}"
           r"\usepackage[active,tightpage]{preview}"
           r"\definecolor{self}{HTML}{B77B00}\definecolor{ld}{HTML}{0072B2}"
           r"\definecolor{ink}{HTML}{16313F}\begin{document}" + "\n")
    tex += "\n".join(r"\begin{preview}\color{ink}$\displaystyle " + f +
                      r"$\end{preview}" for f in FORMULAS)
    tex += "\n" + r"\end{document}" + "\n"
    (OUT / "equations.tex").write_text(tex)
    subprocess.run([TECTONIC] + (["--only-cached"] if args.offline else []) +
                   ["--outdir", str(tmp), str(OUT / "equations.tex")], check=True)

    ns = "http://www.w3.org/2000/svg"
    main = ET.parse(str(tmp / "base.svg"))
    root = main.getroot()
    root.set("width", "1600")
    root.set("height", "1000")
    root.set("viewBox", "0 0 1600 1000")
    root.set("role", "img")
    title = ET.Element(f"{{{ns}}}title")
    title.text = "Figure 28a. LD score and a thresholded count of LD friends"
    desc = ET.Element(f"{{{ns}}}desc")
    desc.text = ("Illustrative, exact population squared genotype correlations among 12 SNPs. "
                 "SNP 4 has squared correlations 0.04, 0.36, 0.64, 1.00, 0.81, 0.49, 0.01, "
                 "and zero with SNPs 8 to 12. Its LD score is 3.35: 1 from itself plus 2.35 "
                 "from other SNPs. At an illustrative r-squared cutoff of 0.20, four other "
                 "SNPs pass. The threshold is not used in computing the LD score. "
                 "GCTA LD friends uses significance of LD, not a universal r-squared cutoff.")
    root.insert(0, desc)
    root.insert(0, title)

    for page, (x, y, width, height) in enumerate([(805, 241, 645, 78), (819, 707, 625, 64)], 1):
        svg_path = tmp / f"equation_{page}.svg"
        subprocess.run([PDFTOCAIRO, "-f", str(page), "-l", str(page), "-svg",
                        str(tmp / "equations.pdf"), str(svg_path)], check=True)
        eq = ET.parse(str(svg_path)).getroot()
        # Cairo glyph identifiers are global within an SVG document.
        ids = {e.get("id"): f"eq{page}_{e.get('id')}" for e in eq.iter() if e.get("id")}
        for e in eq.iter():
            for key, value in list(e.attrib.items()):
                if key == "id":
                    e.set(key, ids[value])
                elif value.startswith("#") and value[1:] in ids:
                    e.set(key, "#" + ids[value[1:]])
                else:
                    e.set(key, re.sub(r"url\(#([^)]*)\)", lambda m: "url(#"+ids.get(m[1],m[1])+")", value))
        view = list(map(float, eq.get("viewBox").split()))
        scale = min(width/view[2], height/view[3])
        g = ET.SubElement(root, f"{{{ns}}}g", transform=f"translate({x},{y}) scale({scale})")
        for child in eq:
            g.append(child)
        number = ET.SubElement(root, f"{{{ns}}}text", x="1494", y=str(y+35),
                               fill="#53656F", **{"font-family":"Calibri, Arial, sans-serif", "font-size":"21", "text-anchor":"end"})
        number.text = (f"({args.equation_number})" if page == 1 else f"({args.equation_number}a)") if args.equation_number else ""

    svg = OUT / (NAME + ".svg")
    main.write(str(svg), encoding="utf-8", xml_declaration=True)
    png = OUT / (NAME + ".png")
    js = "const sharp=require(process.argv[1]);sharp(process.argv[2],{density:144}).resize(3200,2000).png().toFile(process.argv[3]).then(()=>{});"
    subprocess.run([NODE, "-e", js, SHARP, str(svg), str(png)], check=True)
    print(svg)
    print(png)
