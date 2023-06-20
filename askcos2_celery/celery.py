import os
from celery import Celery

# Define readable names for celery workers for status reporting
READABLE_NAMES = {
    "generic_worker": "Generic Worker"
}

# Note: cannot use guest for authenticating with broker unless on localhost
redis_host = os.environ.get("REDIS_HOST", "0.0.0.0")
redis_port = os.environ.get("REDIS_PORT", "6379")
redis_password = os.environ.get("REDIS_PASSWORD", "")
if redis_password:
    redis_password = ":{0}@".format(redis_password)
redis_url = "redis://{password}{host}:{port}".format(
    password=redis_password,
    host=redis_host,
    port=redis_port,
)

rabbit_host = os.environ.get("RABBITMQ_HOST", "0.0.0.0")
rabbit_port = os.environ.get("RABBITMQ_PORT", "5672")
rabbit_username = os.environ.get("RABBITMQ_USERNAME", "guest")
rabbit_password = os.environ.get("RABBITMQ_PASSWORD", "guest")
rabbit_url = "amqp://{username}:{password}@{host}:{port}/".format(
    username=rabbit_username,
    password=rabbit_password,
    host=rabbit_host,
    port=rabbit_port,
)

app = Celery(
    main="askcos2_celery",
    broker=rabbit_url,
    backend=redis_url,
)

# Using a string here means the worker don't have to serialize
# the configuration object to child processes.
app.config_from_object("askcos2_celery.celery_config")


if __name__ == "__main__":
    app.start()
