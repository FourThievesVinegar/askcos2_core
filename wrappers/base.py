import requests
from celery.result import AsyncResult
from pydantic import BaseModel


class BaseWrapper:
    DEFAULT_PREDICTION_URL: str = None
    input_class: type[BaseModel]
    output_class: type[BaseModel]

    name: str = "base"
    prefixes: list[str] = None
    methods_to_bind: dict[str, list[str]] = {
        "get_config": ["GET"],
        "get_doc": ["GET"],
        "call_sync": ["POST"],
        "call_async": ["POST"],
        "retrieve": ["GET"]
    }

    def __init__(self, config: dict):
        self.config = config
        self.prediction_url = config.get(
            "custom_prediction_url",
            self.DEFAULT_PREDICTION_URL
        )
        self.config["prediction_url"] = self.prediction_url
        self.session_sync = requests.Session()

    def get_config(self) -> dict:
        return self.config

    def get_doc(self) -> str:
        return self.__doc__

    def call_sync(self, input: BaseModel) -> BaseModel:
        response = self.session_sync.post(
            self.prediction_url,
            json=input.dict(),
            timeout=self.config["timeout"]
        )
        output = response.json()
        output = self.output_class(**output)

        return output

    async def call_async(self, input: BaseModel, priority: int = 0) -> str:
        from askcos2_celery.tasks import base_task
        async_result = base_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> BaseModel | None:
        result = AsyncResult(task_id)
        output = result.result
        output = self.output_class(**output)

        return output
