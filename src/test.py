import pytest

from main import app
from fastapi.testclient import TestClient

from database import Base, engine


@pytest.fixture()
def test_db():
    Base.metadata.create_all(bind=engine)
    yield
    Base.metadata.drop_all(bind=engine)


client = TestClient(app)


def test_empty_molecule(test_db):
    response = client.get("/substructure")
    assert response.status_code == 404


def test_empty_table(test_db):
    response = client.get("/substructure/CCO")
    assert response.status_code == 200
    assert response.json() == []


def test_one_true(test_db):
    client.post("/molecules", json={"id": 333, "smiles": "Cc1ccccc1"})
    response = client.get("/substructure/c1ccccc1")
    assert response.status_code == 200
    assert response.json() == ["Cc1ccccc1"]


def test_one_false(test_db):
    client.post("/molecules", json={"id": 333, "smiles": "Cc1ccccc1"})
    response = client.get("/substructure/C(=O)O")
    assert response.status_code == 200
    assert response.json() == []


def test_many_true(test_db):
    client.post("/molecules", json={"id": 333, "smiles": "Cc1ccccc1"})
    client.post("/molecules", json={"id": 334, "smiles": "CCO"})
    response = client.get("/substructure/C")
    assert response.status_code == 200
    assert response.json() == ["Cc1ccccc1", "CCO"]


def test_many_false(test_db):
    client.post("/molecules", json={"id": 333, "smiles": "Cc1ccccc1"})
    client.post("/molecules", json={"id": 334, "smiles": "CCO"})
    response = client.get("/substructure/CC(=O)Oc1ccccc1C(=O)O")
    assert response.status_code == 200
    assert response.json() == []


def test_many(test_db):
    client.post("/molecules", json={"id": 333, "smiles": "Cc1ccccc1"})
    client.post("/molecules", json={"id": 334, "smiles": "CCO"})
    response = client.get("/substructure/c1ccccc1")
    assert response.status_code == 200
    assert response.json() == ["Cc1ccccc1"]


def test_invalid_smiles(test_db):
    response = client.get("/substructure/bla")
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid substructure"}
