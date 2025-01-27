from typing import Optional

from pydantic import BaseModel


class MoleculeSchema(BaseModel):
    id: int
    smiles: str


class UpdateMolecule(BaseModel):
    smiles: str


class TaskResultBase(BaseModel):
    task_id: str
    status: str
    result: Optional[str] = None

    class Config:
        from_attributes = True
