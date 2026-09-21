# Complete Lyapunov Functions and Chain-Recurrent Partitions for a Cancer-State Ordinary Differential Equation

**Thesis #17.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

Can a complete Lyapunov construction (mesh-free / radial-basis ideas in the spirit of Argáez, Giesl, and Hafstein) partition a low-dimensional cancer-state ODE into chain-recurrent versus transient regions in a way that is reproducible from the vector field alone, without interpreting basins as treatment response?

The calculation uses one planar tumour–effector field with four declared constants. The coexistence state is an unstable focus, surrounded by a periodic orbit of period 19.88. A Wendland collocation asks the orbital derivative along a speed-normalised copy of the field to equal −1, then rebuilds once with target 0 on the nodes whose flow-aligned stencil mean exceeds the selected cut γ = −0.55.

On a 34×34 grid, 339 of 1156 nodes fail that cut, and the periodic orbit lies within two spacings of a failing node. Along the orbit the rebuilt function changes by 0.024 from peak to peak. A transient from (0.55, 0.88) drops by 1.31 over ten periods and ends in that range. The same orbit coverage appears on a 22×22 grid and on a 42×42 grid. The unstable focus is marked on the finer grids and missed on the coarse one. The failing set is a collocation defect, not a certified chain-recurrent set: the orbital derivative along the orbit is still positive on 46 percent of samples. The window has one attractor and the vector field has no treatment parameter, so the partition is not a response class.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/lyapunov.py` | Seed-labelled collocation (run label 20260921; the arithmetic is deterministic) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Phase portrait, orbital derivative, partition, level sets, decrease, return |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/lyapunov.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the collocation. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Complete Lyapunov functions and chain-recurrent partitions for a cancer-state ordinary differential equation [Internet]. Thesis #17 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-17-complete-lyapunov-cancer-ode

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
