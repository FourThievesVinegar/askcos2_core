import requests
from celery import shared_task
from celery.result import AsyncResult
from pydantic import BaseModel


class BaseWrapper:
    DEFAULT_PREDICTION_URL: str = None
    DEFAULT_MANAGEMENT_URL: str = None
    input_class: type[BaseModel]
    output_class: type[BaseModel]
    hparams_class: type[BaseModel]

    name: str = "base"
    prefixes: list[str] = None
    methods_to_bind: list[str] = [
        "get_hparams", "set_hparams", "call_sync", "call_async", "retrieve"
    ]

    def __init__(self, config: dict, prediction_url: str = None,
                 management_url: str = None):
        self.config = config
        self.prediction_url = prediction_url or self.DEFAULT_PREDICTION_URL
        self.management_url = management_url or self.DEFAULT_MANAGEMENT_URL
        self.session_sync = requests.Session()

    def get_hparams(self) -> BaseModel:
        response = self.session_sync.get(f"{self.management_url}/get_hparams")
        hparams = response.json()
        hparams = self.hparams_class(**hparams)

        return hparams

    def set_hparams(self, hparams: BaseModel) -> bool:
        response = self.session_sync.post(
            f"{self.management_url}/set_hparams",
            json=hparams.dict()
        )

        return response.status_code == 200

    def call_sync(self, input: BaseModel) -> BaseModel:
        response = self.session_sync.post(
            self.prediction_url,
            json=input.dict(),
            timeout=self.config["timeout"]
        )
        output = response.json()
        output = self.output_class(**output)

        return output

    def call_async(self, input: BaseModel, priority: int
                   ) -> tuple[str, BaseModel]:
        async_result = base_task.apply_async(
            self, input=input, priority=priority)
        task_id = async_result.id

        return task_id, input

    @staticmethod
    def retrieve(task_id: str) -> BaseModel | None:
        result = AsyncResult(task_id)
        output = result.result

        return output


class TorchServeWrapper(BaseWrapper):
    name: str = "torchserve"
    methods_to_bind = ["call_sync", "call_async", "retrieve"]

    def get_hparams(self):
        raise NotImplementedError

    def set_hparams(self, hparams: BaseModel):
        raise NotImplementedError


@shared_task
def base_task(wrapper: BaseWrapper, input: BaseModel) -> BaseModel:
    return wrapper.call_sync(input)
