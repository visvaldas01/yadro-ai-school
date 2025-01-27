import json
from os import getenv

import redis
from fastapi import FastAPI, HTTPException, Depends
from rdkit import Chem

from models import Molecule, MoleculeEncoder
from database import SessionLocal, Base, engine
from schemas import UpdateMolecule, MoleculeSchema

import logging

Base.metadata.create_all(bind=engine)

logger = logging.getLogger('uvicorn')

app = FastAPI()

redis_client = redis.Redis(host='localhost', port=6379, db=0)


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def get_cached_result(key: str):
    result = redis_client.get(key)
    if result:
        return json.loads(result)
    return None


def set_cache(key: str, value: dict, expiration: int = 60):
    redis_client.setex(key, expiration, json.dumps(value, cls=MoleculeEncoder))


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/molecules")
def add_molecule(molecule: MoleculeSchema, db: SessionLocal = Depends(get_db)):
    check = db.get(Molecule, molecule.id)
    if check:
        logger.error(f"add_molecule: Molecule {check} already exists")
        raise HTTPException(status_code=400, detail="Molecule already exists")
    rd_molecule = Chem.MolFromSmiles(molecule.smiles)
    if not rd_molecule:
        logger.error(f"add_molecule: {molecule.smiles} is invalid SMILES")
        raise HTTPException(status_code=400, detail="Invalid molecule")
    new_molecule = Molecule(id=molecule.id, smiles=molecule.smiles)
    db.add(new_molecule)
    db.commit()
    logger.info(f"add_molecule: Molecule {new_molecule} added successfully")
    return "Molecule added successfully"


@app.get("/molecules/{molecule_id}")
def get_molecule(molecule_id: int, db: SessionLocal = Depends(get_db)):
    molecule = db.get(Molecule, molecule_id)
    if not molecule:
        logger.error(f"get_molecule: Molecule with id {molecule_id} not found")
        raise HTTPException(status_code=404, detail="Molecule not found")
    logger.info(f"get_molecule: {molecule}")
    return molecule


@app.get("/molecules")
def get_molecules(offset: int = 0, limit: int = 10, db: SessionLocal = Depends(get_db)):
    def molecule_iterator(offset, limit):
        query = db.query(Molecule).offset(offset).limit(limit)
        for molecule in query.all():
            yield molecule

    molecules = list(molecule_iterator(offset, limit))
    return molecules


@app.get("/substructure/{substr}")
def substructure_search(substr: str, db: SessionLocal = Depends(get_db)):
    cache_key = f"search:{substr}"
    cached_result = get_cached_result(cache_key)

    if cached_result:
        response = {"source": "cache", "data": cached_result}
        logger.info(f"substructure_search: {response}")
        return response
    else:
        rd_substr = Chem.MolFromSmiles(substr)
        if rd_substr is None:
            logger.error(f"substructure_search: {substr} is invalid SMILES")
            raise HTTPException(status_code=400, detail="Invalid substructure")
        molecules = [molecule for molecule in db.query(Molecule).all()]
        found_substr = [molecule for molecule in molecules
                        if Chem.MolFromSmiles(molecule.smiles).HasSubstructMatch(rd_substr)]

        search_result = {"query": cache_key, "result": found_substr}
        set_cache(cache_key, search_result)
        response = {"source": "database", "data": search_result}
        logger.info(f"substructure_search: {response}")
        return response


@app.delete("/molecules/{molecule_id}")
def delete_molecule(molecule_id: int, db: SessionLocal = Depends(get_db)):
    check = db.get(Molecule, molecule_id)
    if not check:
        logger.error(f"delete_molecule: Molecule with id {molecule_id} not found")
        raise HTTPException(status_code=404, detail="Molecule not found")
    db.delete(check)
    db.commit()
    logger.info(f"delete_molecule: Molecule {check} deleted successfully")
    return "Molecule deleted successfully"


@app.put("/molecules/{molecule_id}")
def update_molecule(molecule_id: int, molecule: UpdateMolecule, db: SessionLocal = Depends(get_db)):
    check = db.get(Molecule, molecule_id)
    if not check:
        logger.error(f"update_molecule: Molecule with id {molecule_id} not found")
        raise HTTPException(status_code=404, detail="Molecule not found")
    rd_molecule = Chem.MolFromSmiles(molecule.smiles)
    if not rd_molecule:
        logger.error(f"update_molecule: {molecule.smiles} is invalid SMILES")
        raise HTTPException(status_code=400, detail="Invalid molecule")
    check.smiles = molecule.smiles
    db.commit()
    logger.info(f"update_molecule: Molecule {molecule} updated successfully")
    return "Molecule updated successfully"
