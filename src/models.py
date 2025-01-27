import json

from sqlalchemy import Column, Integer, String, Text
from database import Base


class Molecule(Base):
    __tablename__ = 'molecules'
    id = Column(Integer, primary_key=True, index=True)
    smiles = Column(String, index=True)

    def __repr__(self):
        molecule_dict = {'id': self.id, 'smiles': self.smiles}
        return molecule_dict.__repr__()


class MoleculeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Molecule):
            return {'id': obj.id, 'smiles': obj.smiles}
        return json.JSONEncoder.default(self, obj)


class TaskResult(Base):
    __tablename__ = "task_results"

    id = Column(Integer, primary_key=True, index=True)
    task_id = Column(String, unique=True, index=True)
    status = Column(String, index=True)
    result = Column(Text, nullable=True)
