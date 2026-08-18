"""
LigandScope — a reproducible structure-based ligand-binding & property-profiling tool.

Demonstrated on a real anti-virulence docking study against Streptococcus pyogenes
(SpeB / SpeA virulence factors), but the profiling engine is general: paste any SMILES
and it computes the same medicinal-chemistry descriptors used to triage compounds
before docking.

- Docking affinities: REAL AutoDock Vina 1.2.5 output (seed 42) from my thesis pipeline.
- Physicochemical descriptors: computed LIVE by RDKit (reproducible from SMILES).
- Nothing is fabricated; limitations are stated on the 'Methods' tab.

Author: Hafiza Muntaha Fatima  ·  github.com/hafizamuntaha-fatima
"""

from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from src.chem import profile_molecule, draw_molecule, profile_to_dict

DATA = Path(__file__).parent / "data"

st.set_page_config(page_title="LigandScope", page_icon="🧬", layout="wide")


# ---------- data loading (cached so it only reads once) ----------
@st.cache_data
def load_data():
    docking = pd.read_csv(DATA / "docking_results.csv")
    ligands = pd.read_csv(DATA / "ligands.csv")
    return docking, ligands


docking, ligands = load_data()

# target metadata (documented from the study — for display/annotation)
TARGETS = {
    "SpeB": {"pdb": "2UZJ", "note": "Cysteine protease; docked at the catalytic "
             "Cys47/His195 region marked by the bound E64 inhibitor."},
    "SpeA": {"pdb": "1B1Z", "note": "Superantigen exotoxin; blind dock over the surface."},
}


# ---------- header ----------
st.title("🧬 LigandScope")
st.caption(
    "Reproducible structure-based ligand-binding & property profiling · "
    "demo: plant anti-virulence compounds vs *Streptococcus pyogenes* SpeB/SpeA"
)

tab_overview, tab_ligand, tab_any, tab_methods = st.tabs(
    ["📊 Docking overview", "🔬 Ligand profiler", "⚗️ Profile any molecule", "📖 Methods & limitations"]
)


# ================= TAB 1: DOCKING OVERVIEW =================
with tab_overview:
    st.subheader("Binding-affinity ranking (real AutoDock Vina output)")
    st.markdown(
        "Lower (more negative) = stronger predicted binding. "
        "Values are the best affinity per ligand–target pair."
    )

    ranked = docking.sort_values("best_affinity_kcal_mol").reset_index(drop=True)
    ranked.index += 1

    col1, col2 = st.columns([1, 1])
    with col1:
        st.dataframe(
            ranked.rename(columns={
                "receptor": "Target", "pdb_id": "PDB", "ligand": "Ligand",
                "best_affinity_kcal_mol": "Best affinity (kcal/mol)",
            })[["Target", "PDB", "Ligand", "Best affinity (kcal/mol)"]],
            use_container_width=True,
        )
    with col2:
        fig = px.bar(
            docking, x="ligand", y="best_affinity_kcal_mol", color="receptor",
            barmode="group", labels={"best_affinity_kcal_mol": "Affinity (kcal/mol)",
                                     "ligand": "Ligand", "receptor": "Target"},
            title="Predicted binding affinity by ligand and target",
        )
        fig.update_yaxes(autorange="reversed")  # stronger binding lower = shown higher
        st.plotly_chart(fig, use_container_width=True)

    st.info(
        "**Read-out:** EGCG (green tea) is the strongest predicted binder to both virulence "
        "factors (best: EGCG–SpeB, −7.574 kcal/mol). Quercetin is intermediate; thymoquinone "
        "alone is weakest — consistent with the flavonoid driving direct target binding."
    )


# ================= TAB 2: LIGAND PROFILER =================
with tab_ligand:
    st.subheader("Per-ligand physicochemical + drug-likeness profile")
    choice = st.selectbox("Choose a compound", ligands["ligand"].tolist())
    row = ligands[ligands["ligand"] == choice].iloc[0]

    c1, c2 = st.columns([1, 1.3])
    with c1:
        img = draw_molecule(row["smiles"])
        if img is not None:
            st.image(img, caption=f"{choice} (PubChem CID {int(row['pubchem_cid'])})")
        st.caption(f"**Source:** {row['source_plant']}  ·  **Class:** {row['class']}")
        st.caption(f"`SMILES:` {row['smiles']}")

    with c2:
        prof = profile_molecule(row["smiles"])
        d = profile_to_dict(prof)
        st.markdown("**Computed descriptors (RDKit, live):**")
        st.table(pd.DataFrame({
            "Property": ["Molecular weight (g/mol)", "logP (lipophilicity)",
                         "H-bond donors", "H-bond acceptors", "TPSA (Å²)",
                         "Rotatable bonds", "Aromatic rings"],
            "Value": [d["mol_weight"], d["logp"], d["h_bond_donors"],
                      d["h_bond_acceptors"], d["tpsa"], d["rotatable_bonds"],
                      d["aromatic_rings"]],
        }))
        lip = "✅ Pass" if d["lipinski_pass"] else "⚠️ Fail"
        veb = "✅ Pass" if d["veber_pass"] else "⚠️ Fail"
        st.markdown(f"**Lipinski rule-of-5:** {lip}  "
                    f"({d['lipinski_violations']} violation(s))　　**Veber:** {veb}")

    st.divider()
    st.markdown("**This compound's docking scores:**")
    sub = docking[docking["ligand"] == choice][["receptor", "pdb_id", "best_affinity_kcal_mol"]]
    for _, r in sub.iterrows():
        t = r["receptor"]
        st.markdown(f"- **{t} ({r['pdb_id']})** — {r['best_affinity_kcal_mol']} kcal/mol. "
                    f"{TARGETS.get(t, {}).get('note', '')}")


# ================= TAB 3: PROFILE ANY MOLECULE (general tool) =================
with tab_any:
    st.subheader("Profile any molecule from its SMILES")
    st.markdown(
        "The engine is general — paste any SMILES (e.g. a ligand from your own project) "
        "to get the same descriptor + drug-likeness read-out. Try a known drug or a "
        "natural product from PubChem."
    )
    smi = st.text_input("SMILES string", value="CC(=O)Oc1ccccc1C(=O)O")  # aspirin
    if smi:
        prof = profile_molecule(smi)
        if prof is None:
            st.error("Could not parse that SMILES — please check it.")
        else:
            cc1, cc2 = st.columns([1, 1.3])
            with cc1:
                img = draw_molecule(smi)
                if img is not None:
                    st.image(img)
            with cc2:
                d = profile_to_dict(prof)
                st.table(pd.DataFrame({
                    "Property": ["Molecular weight", "logP", "H-bond donors",
                                 "H-bond acceptors", "TPSA (Å²)", "Rotatable bonds",
                                 "Aromatic rings"],
                    "Value": [d["mol_weight"], d["logp"], d["h_bond_donors"],
                              d["h_bond_acceptors"], d["tpsa"], d["rotatable_bonds"],
                              d["aromatic_rings"]],
                }))
                lip = "✅ Pass" if d["lipinski_pass"] else "⚠️ Fail"
                veb = "✅ Pass" if d["veber_pass"] else "⚠️ Fail"
                st.markdown(f"**Lipinski:** {lip} ({d['lipinski_violations']} violation(s))　"
                            f"**Veber:** {veb}")


# ================= TAB 4: METHODS & LIMITATIONS =================
with tab_methods:
    st.subheader("Methods & honest limitations")
    st.markdown("""
**Docking.** AutoDock Vina 1.2.5, fixed seed 42. Receptors: SpeB (PDB **2UZJ**, box centred on
the bound E64 inhibitor at the Cys47/His195 catalytic site) and SpeA (PDB **1B1Z**, blind dock).
Ligands prepared with Open Babel; values are the best affinity per ligand–target pair.

**Descriptors.** Computed live by **RDKit** from each compound's SMILES (verifiable on PubChem
via the CIDs shown). Reproducible: same SMILES → same numbers.

**Limitations (stated up front):**
- Docking predicts **binding pose and relative affinity, not efficacy**.
- **No synergy is claimed.** The flavonoid + thymoquinone combination is a formulation to be
  tested by **in-vitro checkerboard / FICI** assays — synergy cannot be derived from docking.
- Absolute affinities depend on search box, protonation and scoring function (all fixed for
  reproducibility).
- Wet-lab validation (MIC, checkerboard/FICI) is **planned**, not yet performed.

**Why this tool.** It mirrors the computational half of a structure-based drug-discovery
workflow — ligand preparation, binding-affinity ranking, and drug-likeness triage — as a
reusable, reproducible app. The engine works on any target/ligand set.

*Built by Hafiza Muntaha Fatima · BS Biochemistry, GCWUF · github.com/hafizamuntaha-fatima*
""")
