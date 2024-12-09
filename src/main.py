from os import getenv

from fastapi import FastAPI, HTTPException, Depends
from rdkit import Chem

from models import Molecule
from database import SessionLocal, Base, engine
from schemas import UpdateMolecule, MoleculeSchema

Base.metadata.create_all(bind=engine)

app = FastAPI()


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/molecules")
def add_molecule(molecule: MoleculeSchema, db: SessionLocal = Depends(get_db)):
    db_molecule = db.get(Molecule, molecule.id)
    if db_molecule:
        raise HTTPException(status_code=400, detail="Molecule already exists")
    rd_molecule = Chem.MolFromSmiles(molecule.smiles)
    if not rd_molecule:
        raise HTTPException(status_code=400, detail="Invalid molecule")
    new_molecule = Molecule(id=molecule.id, smiles=molecule.smiles)
    db.add(new_molecule)
    db.commit()
    return "Molecule added successfully"


@app.get("/molecules/{molecule_id}")
def get_molecule(molecule_id: int, db: SessionLocal = Depends(get_db)):
    molecule = db.get(Molecule, molecule_id)
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return molecule


@app.get("/molecules")
def get_molecules(db: SessionLocal = Depends(get_db)):
    molecules = db.query(Molecule).all()
    return molecules


@app.get("/substructure/{substr}")
def substructure_search(substr: str, db: SessionLocal = Depends(get_db)):
    rd_substr = Chem.MolFromSmiles(substr)
    if rd_substr is None:
        raise HTTPException(status_code=400, detail="Invalid substructure")
    molecules = [molecule.smiles for molecule in db.query(Molecule).all()]
    found_substr = [molecule for molecule in molecules
                    if Chem.MolFromSmiles(molecule).HasSubstructMatch(rd_substr)]
    return found_substr


@app.delete("/molecules/{molecule_id}")
def delete_molecule(molecule_id: int, db: SessionLocal = Depends(get_db)):
    db_molecule = db.get(Molecule, molecule_id)
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    db.delete(db_molecule)
    db.commit()
    return "Molecule deleted successfully"


@app.put("/molecules/{molecule_id}")
def update_molecule(molecule_id: int, molecule: UpdateMolecule, db: SessionLocal = Depends(get_db)):
    db_molecule = db.get(Molecule, molecule_id)
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    rd_molecule = Chem.MolFromSmiles(molecule.smiles)
    if not rd_molecule:
        raise HTTPException(status_code=400, detail="Invalid molecule")
    db_molecule.smiles = molecule.smiles
    db.commit()
    return "Molecule updated successfully"
