from os import getenv

from fastapi import HTTPException, Depends
from rdkit import Chem

from database import Base, engine, SessionLocal
from models import Molecule, TaskResult
from schemas import UpdateMolecule, MoleculeSchema, TaskResultBase
from utils import logger, app, get_db, get_cached_result, redis_client
from tasks import search_task

Base.metadata.create_all(bind=engine)


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
    logger.info(f"get_molecules: {molecules}")
    return molecules


@app.get("/substructure/{substr}")
def substructure_search(substr: str):
    cache_key = f"search:{substr}"
    cached_result = get_cached_result(cache_key)

    if cached_result:
        response = {"source": "cache", "data": cached_result}
        logger.info(f"substructure_search: {response}")
        return response
    else:
        task_id = redis_client.get(f"task:{cache_key}")
        if not task_id:
            task = search_task.delay(substr, cache_key)
            task_id = task.id
            redis_client.set(f"task:{cache_key}", task_id)
            return {"task_id": task_id, "status": "processing",
                    "message": "Search started. Please check back later."}
        else:
            return {"task_id": task_id, "status": "processing",
                    "message": "Search is in progress. Please check back later."}


@app.get("/result/{task_id}", response_model=TaskResultBase)
def get_result(task_id: str, db: SessionLocal = Depends(get_db)):
    task_result = db.query(TaskResult).filter(TaskResult.task_id == task_id).first()
    if not task_result:
        return {"task_id": task_id, "status": "PENDING", "result": None}
    return task_result


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
