import os

from celery import Celery

celery = Celery(
    'tasks',
    broker='redis://' + os.getenv('REDIS_HOST') + ':6379/0',
    backend='redis://' + os.getenv('REDIS_HOST') + ':6379/0',
    include=['tasks']
)

celery.conf.update(task_track_started=True)
