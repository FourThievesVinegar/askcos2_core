from configs import db_config
from fastapi import Depends, HTTPException, Response, status
from jose import JWTError, jwt
from pydantic import BaseModel
from pymongo import errors, MongoClient
from typing import Annotated, Any
from utils import register_util
from utils.oauth2 import (
    ALGORITHM,
    credentials_exception,
    oauth2_scheme,
    pwd_context,
    SECRET_KEY
)


class User(BaseModel):
    username: str
    hashed_password: str
    email: str | None = None
    full_name: str | None = None
    disabled: bool = False
    is_superuser: bool = False


@register_util(name="user_controller")
class UserController:
    prefixes = ["user"]
    methods_to_bind: dict[str, list[str]] = {
        "register": ["POST"],
        "update": ["POST"],
        "delete": ["DELETE"],
        "enable": ["GET"],
        "disable": ["GET"],
        "get_all_users": ["GET"]
    }

    def __init__(self, util_config: dict[str, Any]):
        """
        Initialize database connection.
        """
        assert util_config["engine"] == "db", \
            "Only the 'db' engine is supported for user controller!"

        self.client = MongoClient(serverSelectionTimeoutMS=1000, **db_config.MONGO)

        database = util_config["database"]
        collection = util_config["collection"]

        try:
            self.client.server_info()
        except errors.ServerSelectionTimeoutError:
            raise ValueError("Cannot connect to mongodb to load user data")
        else:
            self.collection = self.client[database][collection]
            self.db = self.client[database]

    def get_user_by_name(self, username: str) -> User:
        query = {"username": username}
        user = self.collection.find_one(query)
        user = User(**user)

        return user

    def get_current_user(self, token: Annotated[str, Depends(oauth2_scheme)]
                         ) -> User:
        try:
            payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
            username: str = payload.get("sub")
            if username is None:
                raise credentials_exception
        except JWTError:
            raise credentials_exception

        user = self.get_user_by_name(username=username)
        if user is None:
            raise credentials_exception
        if user.disabled:
            raise HTTPException(status_code=400, detail="Disabled user")

        return user

    def register(
        self,
        username: str,
        password: str,
        email: str = None,
        full_name: str = None,
        disabled: bool = False,
        is_superuser: bool = False
    ) -> Response:
        if self.collection.find_one({"username": username}):
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail=f"username: {username} already exists!",
                headers={"WWW-Authenticate": "Bearer"},
            )

        hashed_password = pwd_context.hash(password)
        doc = {
            "username": username,
            "hashed_password": hashed_password,
            "email": email,
            "full_name": full_name,
            "disabled": disabled,
            "is_superuser": is_superuser
        }
        User(**doc)         # data validation, if any
        self.collection.insert_one(doc)

        return Response(content=f"Successfully register user: {username}!")

    def update(
        self,
        username: str,
        password: str,
        email: str = None,
        full_name: str = None,
        disabled: bool = False,
        is_superuser: bool = False,
        token: Annotated[str, Depends(oauth2_scheme)] = None
    ) -> Response:
        user = self.get_current_user(token)
        if not user.is_superuser:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="update operation only permitted by superusers",
                headers={"WWW-Authenticate": "Bearer"},
            )
        hashed_password = pwd_context.hash(password)
        self.collection.update_one(
            {"username": username},
            {"$set": {
                "hashed_password": hashed_password,
                "email": email,
                "full_name": full_name,
                "disabled": disabled,
                "is_superuser": is_superuser
            }}
        )

        return Response(content=f"Successfully update records for {username}!")

    def delete(
        self,
        username: str,
        token: Annotated[str, Depends(oauth2_scheme)]
    ) -> Response:
        user = self.get_current_user(token)
        if not user.is_superuser:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="delete operation only permitted by superusers",
                headers={"WWW-Authenticate": "Bearer"},
            )
        self.collection.delete_one({"username": username})

        return Response(content=f"Successfully delete user: {username}!")

    def enable(
        self,
        username: str,
        token: Annotated[str, Depends(oauth2_scheme)]
    ) -> Response:
        user = self.get_current_user(token)
        if not user.is_superuser:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="enable operation only permitted by superusers",
                headers={"WWW-Authenticate": "Bearer"},
            )
        self.collection.update_one(
            {"username": username},
            {"$set": {"disabled": False}}
        )

        return Response(content=f"Successfully enable user: {username}!")

    def disable(
        self,
        username: str,
        token: Annotated[str, Depends(oauth2_scheme)]
    ) -> Response:
        user = self.get_current_user(token)
        if not user.is_superuser:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="disable operation only permitted by superusers",
                headers={"WWW-Authenticate": "Bearer"},
            )
        self.collection.update_one(
            {"username": username},
            {"$set": {"disabled": True}}
        )

        return Response(content=f"Successfully disable user: {username}!")

    def get_all_users(self, token: Annotated[str, Depends(oauth2_scheme)]
                      ) -> list[User]:
        user = self.get_current_user(token)
        if not user.is_superuser:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="get_all_users operation only permitted by superusers",
                headers={"WWW-Authenticate": "Bearer"},
            )
        users = self.collection.find()
        users = [User(**u) for u in users]

        return users
