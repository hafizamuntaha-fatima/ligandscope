"""
chem.py — molecular property + drug-likeness engine for LigandScope.

Every value here is computed live by RDKit from the molecule's SMILES string.
Nothing is hard-coded, so the numbers are reproducible: give the same SMILES,
get the same descriptors. These are standard medicinal-chemistry descriptors
(the same ones used to triage compounds before docking).

Author: Hafiza Muntaha Fatima
"""

from __future__ import annotations
from dataclasses import dataclass, asdict

from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, Crippen, Lipinski, rdMolDescriptors


@dataclass
class MolProfile:
    """A compound's computed physicochemical profile."""
    smiles: str
    mol_weight: float          # g/mol
    logp: float                # Crippen logP (lipophilicity)
    h_bond_donors: int         # number of H-bond donors
    h_bond_acceptors: int      # number of H-bond acceptors
    tpsa: float                # topological polar surface area (Å^2)
    rotatable_bonds: int       # molecular flexibility
    aromatic_rings: int

    # rule-based drug-likeness verdicts (True = passes / drug-like)
    lipinski_pass: bool        # Lipinski rule-of-5 (<=1 violation)
    lipinski_violations: int
    veber_pass: bool           # Veber: rot. bonds <=10 and TPSA <=140


def mol_from_smiles(smiles: str):
    """Parse a SMILES string into an RDKit molecule (or None if invalid)."""
    return Chem.MolFromSmiles(smiles)


def profile_molecule(smiles: str) -> MolProfile | None:
    """
    Compute the physicochemical + drug-likeness profile of one molecule.

    Lipinski rule-of-5 (oral drug-likeness): MW <= 500, logP <= 5,
    H-bond donors <= 5, H-bond acceptors <= 10. Up to one violation is
    still considered 'drug-like'.
    Veber rules: rotatable bonds <= 10 AND TPSA <= 140 Å^2.
    """
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rot = Lipinski.NumRotatableBonds(mol)
    ar_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

    # count Lipinski violations
    violations = 0
    if mw > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if hbd > 5:
        violations += 1
    if hba > 10:
        violations += 1

    return MolProfile(
        smiles=smiles,
        mol_weight=round(mw, 2),
        logp=round(logp, 2),
        h_bond_donors=hbd,
        h_bond_acceptors=hba,
        tpsa=round(tpsa, 2),
        rotatable_bonds=rot,
        aromatic_rings=ar_rings,
        lipinski_pass=(violations <= 1),
        lipinski_violations=violations,
        veber_pass=(rot <= 10 and tpsa <= 140),
    )


def draw_molecule(smiles: str, size: tuple[int, int] = (350, 300)):
    """Return a 2D structure image (PIL) for the given SMILES, or None."""
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    return Draw.MolToImage(mol, size=size)


def profile_to_dict(profile: MolProfile) -> dict:
    """Flatten a MolProfile to a plain dict (for tables/JSON)."""
    return asdict(profile)
