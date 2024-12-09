from pydantic import BaseModel


class MoleculeSchema(BaseModel):
    id: int
    smiles: str


class UpdateMolecule(BaseModel):
    smiles: str
