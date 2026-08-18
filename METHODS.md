# LigandScope — Methods Note (1 page)

**Author:** Hafiza Muntaha Fatima · BS Biochemistry, Government College Women University Faisalabad

## Purpose
LigandScope packages the **computational half of a structure-based drug-discovery workflow** —
ligand preparation, binding-affinity ranking, and physicochemical/drug-likeness triage — into a
single reproducible, interactive tool. It is demonstrated on a real anti-virulence study but the
engine is target-agnostic.

## Study system (demonstration)
A plant-based **anti-virulence** hypothesis against *Streptococcus pyogenes* (Group A
*Streptococcus*): rather than killing the bacterium, disable its secreted virulence factors.
- **Targets:** SpeB (cysteine protease, PDB **2UZJ**) and SpeA (superantigen exotoxin, PDB **1B1Z**).
- **Ligands:** EGCG (green tea), quercetin (flavonol), thymoquinone (*Nigella sativa*).

## Computational methods
1. **Receptor preparation** — PDB structures cleaned; for SpeB the search box was centred on the
   bound **E64** inhibitor, i.e. the catalytic **Cys47/His195** site; SpeA was blind-docked.
2. **Ligand preparation** — 3D structures from SMILES (PubChem CIDs), Open Babel.
3. **Docking** — **AutoDock Vina 1.2.5**, fixed **seed 42** for reproducibility; best affinity
   (kcal/mol) recorded per ligand–target pair.
4. **Property profiling** — **RDKit** computes MW, logP, H-bond donors/acceptors, TPSA, rotatable
   bonds; **Lipinski rule-of-5** and **Veber** rules applied as drug-likeness filters.

## Key result
| Rank | Target | Ligand | Best affinity (kcal/mol) |
|---|---|---|---|
| 1 | SpeB (2UZJ) | **EGCG** | **−7.574** |
| 2 | SpeA (1B1Z) | EGCG | −7.497 |
| 3 | SpeA (1B1Z) | Quercetin | −7.123 |

EGCG is the strongest predicted binder to both virulence factors; the flavonoid scaffold appears
to drive direct target engagement.

## Limitations (stated explicitly)
- Docking gives **pose + relative affinity, not efficacy**.
- **No synergy is claimed**; the flavonoid + thymoquinone combination is a formulation for
  **in-vitro checkerboard / FICI** testing.
- Wet-lab validation (MIC, checkerboard/FICI) is **planned**, not performed.

## Relevance to structure-based drug design (MD3 lab context)
This is the same computational entry point used before experimental structure determination and
affinity optimisation: rank candidate ligands, filter for drug-likeness, then prioritise for
biophysical validation. The tool can be pointed at any protein target and ligand series.

*Reproduce: `pip install -r requirements.txt && streamlit run app.py`. Docking regenerates from
`notebooks/molecular_docking_pipeline.ipynb` in the spyogenes-docking repo.*
