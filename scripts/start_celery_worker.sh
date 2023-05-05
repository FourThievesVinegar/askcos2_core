celery -A askcos2_celery worker -c "$1" -Q "$2" -n "$3"@%h --pool=gevent
