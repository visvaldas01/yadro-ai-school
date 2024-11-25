from pydantic import BaseModel


class Molecule(BaseModel):
    id: str
    smiles: str

class UpdateMolecule(BaseModel):
    smiles: str