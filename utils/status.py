from configs import db_config
from typing import Any
from utils import register_util
from pydantic import BaseModel
from pymongo import errors, MongoClient


class CollectionData(BaseModel):
    name: str
    description: str
    url: str
    total: int
    details: dict[str, int | None]
    field: str


class StatusDatabaseResponse(BaseModel):
    collections: list[CollectionData]


@register_util(name="status_database")
class StatusDatabase:
    """
    Util class for showing database status
    """
    prefixes = ["status/database"]
    methods_to_bind: dict[str, list[str]] = {
        "get": ["GET"]
    }

    def __init__(self, util_config: dict[str, Any] = None):
        self.client = MongoClient(serverSelectionTimeoutMS=1000, **db_config.MONGO)
        database = "askcos"

        try:
            self.client.server_info()
        except errors.ServerSelectionTimeoutError:
            raise ValueError("Cannot connect to mongodb for api logging")
        else:
            self.db = self.client[database]

        self.collections = [
            {
                "name": "buyables",
                "description": "Collection of prices for available building blocks.",
                "url": "/buyables/",
                "total": 0,
                "details": {},
                "field": "Source"
            },
            {
                "name": "chemicals",
                "description": "Collection of chemical popularity in training data.",
                "url": "",
                "total": 0,
                "details": {},
                "field": "Template Set"
            },
            {
                "name": "reactions",
                "description":
                    "Collection of training reaction data for retro templates.",
                "url": "",
                "total": 0,
                "details": {},
                "field": "Template Set"
            },
            {
                "name": "retro_templates",
                "description": "Collection of retrosynthetic templates.",
                "url": "",
                "total": 0,
                "details": {},
                "field": "Template Set"
            }
        ]

    def get(self) -> StatusDatabaseResponse:
        """Handle GET requests for database status."""
        data = []
        for collection in self.collections:
            total, details = self._get_counts(collection=collection)
            collection["total"] = total
            collection["details"] = details
            data.append(CollectionData(**collection))
        resp = StatusDatabaseResponse(collections=data)

        return resp

    def _get_counts(self, collection: dict[str, Any]
                    ) -> tuple[int, dict[str, int | None]]:
        _field_map = {
            "Source": "$source",
            "Template Set": "$template_set"
        }

        _collection = self.db[collection["name"]]
        estimated_count = _collection.estimated_document_count()
        if estimated_count == collection["total"]:
            # nothing new, return already computed counts
            return collection["total"], collection["details"]
        else:
            # otherwise, recompute
            cursor = _collection.aggregate(
                [
                    {
                        "$group": {
                            "_id": _field_map[collection["field"]],
                            "count": {"$sum": 1},
                        }
                    },
                ]
            )
            details = {x["_id"]: x["count"] for x in cursor}
            total = sum(details.values())

            return total, details
