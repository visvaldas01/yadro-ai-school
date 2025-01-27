from celery.signals import task_postrun
from rdkit import Chem

from celery_worker import celery
from database import SessionLocal
from models import TaskResult, Molecule
from utils import set_cache, redis_client, logger, get_cached_result


@celery.task
def search_task(substr, cache_key):
    session = SessionLocal()
    rd_substr = Chem.MolFromSmiles(substr)
    if rd_substr is None:
        logger.error(f"substructure_search: {substr} is invalid SMILES")
        raise Exception("Invalid substructure")
    molecules = [molecule for molecule in session.query(Molecule).all()]
    session.close()
    found_substr = [molecule for molecule in molecules
                    if Chem.MolFromSmiles(molecule.smiles).HasSubstructMatch(rd_substr)]

    search_result = {"query": cache_key, "result": found_substr}
    set_cache(cache_key, search_result)
    redis_client.delete(f"task:{cache_key}")
    response = {"source": "database", "data": get_cached_result(cache_key)}
    logger.info(f"substructure_search: {response}")
    return response


@task_postrun.connect
def task_postrun_handler(sender=None, task_id=None, task=None, args=None, kwargs=None, retval=None, state=None,
                         **extras):
    session = SessionLocal()
    try:
        task_result = session.query(TaskResult).filter(TaskResult.task_id == task_id).first()
        if not task_result:
            task_result = TaskResult(task_id=task_id)
            session.add(task_result)
        task_result.status = state
        if state == 'SUCCESS':
            task_result.result = str(retval)
        else:
            task_result.result = None
        session.commit()
    except Exception as e:
        session.rollback()
        print(f"Ошибка при сохранении результата задачи: {e}")
    finally:
        session.close()
