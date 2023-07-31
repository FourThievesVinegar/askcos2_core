import requests
from celery.result import AsyncResult
from pydantic import BaseModel
from typing import Literal


class BaseResponse(BaseModel):
    status_code: int
    message: str
    result_format: Literal["json", "base64"]
    result: str | list | dict


class BaseWrapper:
    input_class: type[BaseModel]
    output_class: type[BaseModel]
    response_class: type[BaseResponse]

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
        self.prediction_url = config["deployment"].get("custom_prediction_url", "")
        if not self.prediction_url:
            self.prediction_url = config["deployment"]["default_prediction_url"]

        self.config["prediction_url_in_use"] = self.prediction_url
        self.session_sync = requests.Session()

    def get_config(self) -> dict:
        return self.config

    def get_doc(self) -> str:
        return self.__doc__

    def call_raw(self, input: BaseModel) -> BaseModel:
        response = self.session_sync.post(
            self.prediction_url,
            json=input.dict(),
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = self.output_class(**output)

        return output

    def call_sync(self, input: BaseModel) -> BaseResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: BaseModel, priority: int = 0) -> str:
        from askcos2_celery.tasks import base_task
        async_result = base_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> BaseResponse | None:
        result = AsyncResult(task_id)
        response = result.result
        response = self.response_class(**response)

        return response

    @staticmethod
    def convert_output_to_response(output: BaseModel) -> BaseResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result_format": "json",
            "result": output.results
        }
        response = BaseResponse(**response)

        return response
