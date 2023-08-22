from pydantic import BaseModel


class BaseResponse(BaseModel):
    status_code: int
    message: str
    result: str | int | float | list | dict
