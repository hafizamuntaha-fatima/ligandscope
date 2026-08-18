# 🧬 LigandScope

**A reproducible structure-based ligand-binding & drug-likeness profiling tool.**

Built by **Hafiza Muntaha Fatima** (BS Biochemistry, GCWUF) · [github.com/hafizamuntaha-fatima](https://github.com/hafizamuntaha-fatima)

LigandScope mirrors the *computational half* of a structure-based drug-discovery workflow —
ligand preparation, **binding-affinity ranking**, and **drug-likeness triage** — as a reusable,
reproducible web app. It is demonstrated here on a real anti-virulence docking study against
*Streptococcus pyogenes*, but the profiling engine is **general** (paste any SMILES).

## What it does
- **Docking overview** — ranks real AutoDock Vina affinities for plant compounds (EGCG, quercetin,
  thymoquinone) against the SpeB (PDB 2UZJ) and SpeA (PDB 1B1Z) virulence factors.
- **Ligand profiler** — per-compound 2D structure + RDKit descriptors (MW, logP, HBD/HBA, TPSA,
  rotatable bonds) + Lipinski/Veber drug-likeness verdicts, alongside its docking scores.
- **Profile any molecule** — the descriptor engine runs on any SMILES you paste.
- **Methods & limitations** — honest, up-front (docking ≠ efficacy; no synergy claimed; wet-lab
  validation is planned).

## Data provenance (everything is real)
- **Docking affinities:** AutoDock Vina 1.2.5, fixed seed 42, from my thesis pipeline
  ([spyogenes-docking](https://github.com/hafizamuntaha-fatima/spyogenes-docking)).
- **Structures:** RCSB PDB 2UZJ (SpeB), 1B1Z (SpeA).
- **Ligand SMILES:** PubChem (CIDs in `data/ligands.csv`) — independently verifiable.
- **Descriptors:** computed live by RDKit; same SMILES → same numbers.

## Run locally
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Deploy (Streamlit Community Cloud)
1. Push this repo to GitHub (can be private).
2. On https://share.streamlit.io → **New app** → select this repo → main file `app.py`.
3. In **Advanced settings**, set Python version to **3.11** (RDKit wheels are available there).
4. Deploy → you get a public URL (shareable even though the repo stays private).

## Project layout
```
app.py                 Streamlit UI (4 tabs)
src/chem.py            RDKit descriptor + drug-likeness engine
data/docking_results.csv   real Vina affinities
data/ligands.csv           ligand names, PubChem CIDs, SMILES, source
requirements.txt
```

## Limitations
Docking predicts binding pose and *relative* affinity, not efficacy. No synergy is claimed —
the flavonoid + thymoquinone combination is a formulation to be tested by in-vitro
checkerboard / FICI assays. Wet-lab validation is planned, not yet performed.
