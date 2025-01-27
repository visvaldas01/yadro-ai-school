import json
import logging
import os

import redis
from fastapi import FastAPI

from models import MoleculeEncoder
from database import SessionLocal

logger = logging.getLogger('uvicorn')

app = FastAPI()

redis_client = redis.Redis(host=os.getenv('REDIS_HOST'), port=6379, db=0)


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
