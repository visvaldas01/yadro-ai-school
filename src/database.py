import os

from sqlalchemy import create_engine
from sqlalchemy.engine import URL
from sqlalchemy.orm import declarative_base, sessionmaker

DATABASE_URL = URL.create(
    drivername=os.getenv('POSTGRES_DRIVER', 'postgresql'),
    username=os.getenv('POSTGRES_USER', 'user'),
    password=os.getenv('POSTGRES_PASSWORD', 'password'),
    host=os.getenv('POSTGRES_HOST'),
    database=os.getenv('POSTGRES_DB', 'db')
)

engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()
