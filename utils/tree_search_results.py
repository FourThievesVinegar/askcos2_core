import json
from configs import db_config
from datetime import datetime
from fastapi import Depends, HTTPException, Response
from pydantic import BaseModel, constr, Field
from pymongo import errors, MongoClient
from typing import Annotated, Any, Literal
from utils import register_util
from utils.oauth2 import oauth2_scheme
from utils.registry import get_util_registry

# TODO: fix the time zone issue


class TreeSearchSavedResults(BaseModel):
    """
    Reimplemented SavedResults model
    """
    user: str
    description: constr(max_length=1000) = None
    created: datetime = Field(default_factory=datetime.now)
    modified: datetime = Field(default_factory=datetime.now)
    dt: constr(max_length=200) = None
    result_id: constr(max_length=64) = None
    result: dict | None = None
    setting: dict | None = None
    tags: str | None = None
    check_date: str | None = None
    result_state: constr(max_length=64) = None
    result_type: Literal["tree_builder", "ipp"] = "tree_builder"
    revision: int = 0

    public: bool = False
    shared_with: list[str] = Field(default_factory=list)


@register_util(name="tree_search_results_controller")
class TreeSearchResultsController:
    """
    Util class for tree search results controller
    """
    prefixes = ["results"]
    methods_to_bind: dict[str, list[str]] = {
        "list": ["GET"],
        "retrieve": ["GET"],
        "create": ["POST"],
        "update": ["PUT"],
        "destroy": ["DELETE"],
        "share": ["GET"],
        "unshare": ["GET"],
        "add": ["POST"],
        "remove": ["DELETE"]
    }

    def __init__(self, util_config: dict[str, Any]):
        """
        Initialize database connection.
        """
        self.client = MongoClient(serverSelectionTimeoutMS=1000, **db_config.MONGO)
        database = util_config["database"]
        collection = util_config["collection"]

        try:
            self.client.server_info()
        except errors.ServerSelectionTimeoutError:
            raise ValueError("Cannot connect to mongodb to load results")
        else:
            self.collection = self.client[database][collection]
            self.db = self.client[database]

    def list(self, token: Annotated[str, Depends(oauth2_scheme)]
             ) -> list[TreeSearchSavedResults]:
        """
        API endpoint for accessing user results list.

        Method: GET

        Returns: list of tree builder results belonging to or shared with the currently
            authenticated user
        """
        user_controller = get_util_registry().get_util(module="user_controller")
        user = user_controller.get_current_user(token)

        query = {
            "$or": [
                {"user": user.username},
                {
                    "public": True,
                    "shared_with": user.username
                }
            ]
        }
        # unsetting the result and setting fields to avoid gigantic returns
        cursor = self.collection.aggregate(
            [
                {"$match": query},
                {"$unset": ["result", "setting"]},
                {"$sort": {"dt": -1}}
            ]
        )
        results = []
        for r in cursor:
            result = TreeSearchSavedResults(**r)
            results.append(result)

        return results

    def retrieve(self, result_id: str, token: Annotated[str, Depends(oauth2_scheme)]
                 ) -> TreeSearchSavedResults:
        """
        API endpoint to retrieve specified result by result_id.

        Method: GET
        """

        user_controller = get_util_registry().get_util(module="user_controller")
        user = user_controller.get_current_user(token)

        query = {
            "result_id": result_id,
            "$or": [
                {"user": user.username},
                {
                    "public": True,
                    "shared_with": user.username
                }
            ]
        }
        result = self.collection.find_one(query)
        if not result:
            raise HTTPException(
                status_code=404,
                detail=f"Result with id {result_id} not found or not viewable!"
            )
        if result["result_state"] in ["completed", "ipp"]:
            result["_id"] = str(result["_id"])
            result = TreeSearchSavedResults(**result)

            return result
        else:
            raise HTTPException(
                status_code=409,
                detail=f"Job not yet complete for id: {result_id}!"
            )

    def create(
        self,
        result: TreeSearchSavedResults,
        token: Annotated[str, Depends(oauth2_scheme)]
    ) -> Response:
        """
        API endpoint to create a result.

        Method: POST
        """
        user_controller = get_util_registry().get_util(module="user_controller")
        user = user_controller.get_current_user(token)

        self.collection.insert_one(result.dict())

        return Response(content=f"Successfully create a result for user: "
                                f"{user.username}!")

    def update(
        self,
        result_id: str,
        updated_result: dict,
        token: Annotated[str, Depends(oauth2_scheme)]
    ) -> Response:
        """
        API endpoint to update a result.

        Method: PUT
        """
        user_controller = get_util_registry().get_util(module="user_controller")
        user = user_controller.get_current_user(token)

        query = {
            "result_id": result_id,
            "revision": updated_result["revision"] - 1,
            "$or": [
                {"user": user.username},
                {
                    "public": True,
                    "shared_with": user.username
                }
            ]
        }
        updated_result = {
            k: v for k, v in updated_result.items() if k in [
                "result", "settings", "description", "tags", "modified", "revision"
            ] and v
        }
        res = self.collection.update_one(
            query,
            {"$set": updated_result}
        )
        if not res.matched_count:
            content = {
                "message": f"Result {result_id} not editable by user {user.username},"
                           f"or result changed since initial access!"
            }
            return Response(
                content=json.dumps(content),
                status_code=404,
                media_type="application/json"
            )
        else:
            return Response(content=f"Successfully update result {result_id}!")

    def destroy(self, result_id: str, token: Annotated[str, Depends(oauth2_scheme)]
                ) -> Response:
        """
        API endpoint to delete specified result by result_id.

        Method: DELETE
        """

        user_controller = get_util_registry().get_util(module="user_controller")
        user = user_controller.get_current_user(token)

        query = {
            "user": user.username,
            "result_id": result_id
        }
        try:
            res = self.collection.delete_one(query)
        except Exception:
            content = {
                "message": f"Could not delete result: {result_id}!"
            }
            return Response(
                content=json.dumps(content),
                status_code=500,
                media_type="application/json"
            )
        if not res.deleted_count:
            content = {
                "message": f"Failed to delete result {result_id}!"
            }
            return Response(
                content=json.dumps(content),
                status_code=500,
                media_type="application/json"
            )
        else:
            return Response(content=f"Successfully delete result: {result_id}!")

    def share(self, result_id: str, token: Annotated[str, Depends(oauth2_scheme)]
              ) -> Response:
        """
        API endpoint to make public specific result by result_id.

        Method: GET
        """

        user_controller = get_util_registry().get_util(module="user_controller")
        user = user_controller.get_current_user(token)

        query = {
            "user": user.username,
            "result_id": result_id
        }
        res = self.collection.update_one(
            query,
            {"$set": {"public": True}}
        )
        if not res.matched_count:
            content = {
                "message": f"Result {result_id} not editable by user {user.username}!"
            }
            return Response(
                content=json.dumps(content),
                status_code=401,
                media_type="application/json"
            )
        else:
            return Response(content=f"Successfully make result public: {result_id}!")

    def unshare(self, result_id: str, token: Annotated[str, Depends(oauth2_scheme)]
                ) -> Response:
        """
        API endpoint to make private specific result by result_id.

        Method: GET
        """

        user_controller = get_util_registry().get_util(module="user_controller")
        user = user_controller.get_current_user(token)

        query = {
            "user": user.username,
            "result_id": result_id
        }
        res = self.collection.update_one(
            query,
            {"$set": {"public": False}}
        )
        if not res.matched_count:
            content = {
                "message": f"Result {result_id} not editable by user {user.username}!"
            }
            return Response(
                content=json.dumps(content),
                status_code=401,
                media_type="application/json"
            )
        else:
            return Response(content=f"Successfully make result private: {result_id}!")

    def add(self, result_id: str, token: Annotated[str, Depends(oauth2_scheme)]
            ) -> Response:
        """
        API endpoint to add specified result to the current user by result_id.

        Method: PUT
        """

        user_controller = get_util_registry().get_util(module="user_controller")
        user = user_controller.get_current_user(token)

        query = {
            "result_id": result_id,
            "public": True
        }
        res = self.collection.update_one(
            query,
            {"$push": {"shared_with": user.username}}
        )
        if not res.matched_count:
            content = {
                "message": f"Result {result_id} not found or is not public!"
            }
            return Response(
                content=json.dumps(content),
                status_code=404,
                media_type="application/json"
            )
        else:
            return Response(content=f"Successfully share result {result_id} for user: "
                                    f"{user.username}!")

    def remove(self, result_id: str, token: Annotated[str, Depends(oauth2_scheme)]
               ) -> Response:
        """
        API endpoint to remove specified result to the current user by result_id.

        Method: DELETE
        """

        user_controller = get_util_registry().get_util(module="user_controller")
        user = user_controller.get_current_user(token)

        query = {
            "result_id": result_id,
            "public": True
        }
        res = self.collection.update_one(
            query,
            {"$pull": {"shared_with": user.username}}
        )
        if not res.matched_count:
            content = {
                "message": f"Result {result_id} not found or is not public!"
            }
            return Response(
                content=json.dumps(content),
                status_code=404,
                media_type="application/json"
            )
        else:
            return Response(content=f"Successfully unshare result {result_id} for "
                                    f"user: {user.username}!")

    def update_result_state(
        self,
        result_id: str,
        state: str,
        token: Annotated[str, Depends(oauth2_scheme)] = None
    ) -> None:
        user_controller = get_util_registry().get_util(module="user_controller")
        user = user_controller.get_current_user(token)

        query = {
            "user": user.username,
            "result_id": result_id
        }
        self.collection.update_one(
            query,
            {"$set": {"result_state": state}}
        )

    def save_results(
        self,
        result_id: str,
        result: BaseModel,
        token: Annotated[str, Depends(oauth2_scheme)] = None
    ) -> None:
        user_controller = get_util_registry().get_util(module="user_controller")
        user = user_controller.get_current_user(token)

        query = {
            "user": user.username,
            "result_id": result_id
        }
        result_doc = result.dict()
        if "result_id" not in result_doc or not result_doc["result_id"]:
            result_doc["result_id"] = result_id

        self.collection.update_one(
            query,
            {"$set": {"result": result_doc}}
        )
