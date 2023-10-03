from config import db_config
from pydantic import BaseModel
from pymongo import errors, MongoClient
from typing import Any
from utils import register_util


class ReactionsInput(BaseModel):
    ids: list[int]
    template_set: str = None


class ReactionsResponse(BaseModel):
    reactions: list


@register_util(name="reactions")
class Reactions:
    """Util class for Reactions"""
    prefixes = ["reactions"]
    methods_to_bind: dict[str, list[str]] = {
        "post": ["POST"]
    }

    def __init__(self, util_config: dict[str, Any]):
        self.client = MongoClient(serverSelectionTimeoutMS=1000, **db_config.MONGO)
        database = "askcos"
        collection = "reactions"

        try:
            self.client.server_info()
        except errors.ServerSelectionTimeoutError:
            raise ValueError("Cannot connect to mongodb for reactions")
        else:
            self.db = self.client[database]
            self.collection = self.client[database][collection]

    def post(self, data: ReactionsInput) -> ReactionsResponse:
        query = {"reaction_id": {"$in": data.ids}}
        if data.template_set:
            # Processing for template subsets which use the same historian data
            query["template_set"] = data.template_set.split(":")[0]

        reactions_by_ids = list(self.collection.find(query))
        resp = ReactionsResponse(reactions=reactions_by_ids)

        return resp
