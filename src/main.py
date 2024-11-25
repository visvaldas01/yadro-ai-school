from os import getenv

from fastapi import FastAPI, HTTPException
from rdkit import Chem

from models import Molecule, UpdateMolecule

app = FastAPI()
molecules = {}


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/molecules")
def add_molecule(molecule: Molecule):
    if molecule.id in molecules:
        raise HTTPException(status_code=400, detail="Molecule already exists")
    rd_molecule = Chem.MolFromSmiles(molecule.smiles)
    if rd_molecule is None:
        raise HTTPException(status_code=400, detail="Invalid molecule")
    molecules[molecule.id] = molecule.smiles
    return "Molecule added successfully"


@app.get("/molecules/{molecule_id}")
def get_molecule(molecule_id: int):
    if molecule_id not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    molecule = molecules[molecule_id]
    return {"id": molecule.id, "smiles": molecule.smiles}


@app.get("/molecules")
def get_molecules():
    return molecules


@app.get("/substructure/{substr}")
def substructure_search(substr: str):
    rd_substr = Chem.MolFromSmiles(substr)
    if rd_substr is None:
        raise HTTPException(status_code=400, detail="Invalid substructure")
    return [molecule for molecule in molecules.values()
            if Chem.MolFromSmiles(molecule).HasSubstructMatch(rd_substr)]


@app.delete("/molecules/{molecule_id}")
def delete_molecule(molecule_id: str):
    if molecule_id not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    del molecules[molecule_id]
    return "Molecule deleted successfully"


@app.put("/molecules/{molecule_id}")
def update_molecule(molecule_id: str, molecule: UpdateMolecule):
    if molecule_id not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    rd_molecule = Chem.MolFromSmiles(molecule.smiles)
    if rd_molecule is None:
        raise HTTPException(status_code=400, detail="Invalid molecule")
    molecules[molecule_id] = molecule.smiles
    return "Molecule updated successfully"
