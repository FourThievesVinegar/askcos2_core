import os
import uuid
from configs import db_config
from fastapi import Depends, HTTPException, Response, status
from jose import JWTError, jwt
from keycloak import KeycloakOpenID
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
        "reset_password": ["POST"],
        "delete": ["DELETE"],
        "enable": ["GET"],
        "disable": ["GET"],
        "get_all_users": ["GET"],
        "get_current_user": ["GET"],
        "promote": ["GET"],
        "demote": ["GET"],
        "am_i_superuser": ["GET"]
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

        # create default admin account if there's none
        users = self.collection.find()
        users = [User(**u) for u in users]
        if not any(u.is_superuser for u in users):
            self.register(
                username="admin",
                password="reallybadpassword"
            )

        server_url = os.environ.get("KEYCLOAK_SERVER_URL", "")
        realm_name = os.environ.get("KEYCLOAK_REALM_NAME", "")
        client_id = os.environ.get("KEYCLOAK_CLIENT_ID", "")
        client_secret_key = os.environ.get("KEYCLOAK_CLIENT_SECRET_KEY", "")

        if server_url:
            self.keycloak_openid = KeycloakOpenID(
                server_url=server_url,
                realm_name=realm_name,
                client_id=client_id,
                client_secret_key=client_secret_key
            )
            self.keycloak_public_key = "-----BEGIN PUBLIC KEY-----\n" + \
                                       self.keycloak_openid.public_key() + \
                                       "\n-----END PUBLIC KEY-----"

    def get_user_by_name(self, username: str) -> User | None:
        query = {"username": username}
        try:
            user = self.collection.find_one(query)
            user = User(**user)
        except:
            user = None

        return user

    def get_current_user(self, token: Annotated[str, Depends(oauth2_scheme)]
                         ) -> User:
        try:
            payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
            username: str = payload.get("sub")
            if username is None:
                raise JWTError
            user = self.get_user_by_name(username=username)

        except JWTError:
            # second try for keycloak
            try:
                options = {
                    "verify_signature": True,
                    "verify_aud": False,
                    "verify_exp": False
                }
                token_info = self.keycloak_openid.decode_token(
                    token=token,
                    key=self.keycloak_public_key,
                    options=options
                )
                username = token_info["preferred_username"] + "_sso"
                user = self.get_user_by_name(username=username)
                # automatically register in mongo if user doesn't exist
                if user is None:
                    self.register(
                        username=username,
                        password=str(uuid.uuid4())
                    )
                    user = self.get_user_by_name(username=username)

            except:
                raise credentials_exception

        if user is None:
            raise credentials_exception

        if user.disabled:
            raise HTTPException(status_code=400, detail="Disabled user")

        return user

    def am_i_superuser(self, token: Annotated[str, Depends(oauth2_scheme)]
                       ) -> bool:
        user = self.get_current_user(token=token)

        return user.is_superuser

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
        email: str = None,
        full_name: str = None,
        disabled: bool = False,
        token: Annotated[str, Depends(oauth2_scheme)] = None
    ) -> Response:
        user = self.get_current_user(token)
        if user.username == username or user.is_superuser:
            self.collection.update_one(
                {"username": username},
                {"$set": {
                    "email": email,
                    "full_name": full_name,
                    "disabled": disabled
                }}
            )

        else:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="update operation only permitted by the owner superusers",
                headers={"WWW-Authenticate": "Bearer"},
            )

        return Response(content=f"Successfully update records for {username}!")

    def reset_password(
        self,
        username: str,
        password: str,
        token: Annotated[str, Depends(oauth2_scheme)] = None
    ) -> Response:
        user = self.get_current_user(token)
        if user.username == username or user.is_superuser:
            hashed_password = pwd_context.hash(password)
            self.collection.update_one(
                {"username": username},
                {"$set": {
                    "hashed_password": hashed_password
                }}
            )
        else:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="password reset only permitted by the owner or superusers",
                headers={"WWW-Authenticate": "Bearer"},
            )

        return Response(content=f"Successfully reset the password for {username}!")

    def promote(
        self,
        username: str,
        token: Annotated[str, Depends(oauth2_scheme)] = None
    ) -> Response:
        user = self.get_current_user(token)
        if user.is_superuser:
            self.collection.update_one(
                {"username": username},
                {"$set": {
                    "is_superuser": True
                }}
            )
        else:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="make superuser only permitted by superusers",
                headers={"WWW-Authenticate": "Bearer"},
            )

        return Response(content=f"Successfully reset the password for {username}!")

    def demote(
        self,
        username: str,
        token: Annotated[str, Depends(oauth2_scheme)] = None
    ) -> Response:
        user = self.get_current_user(token)
        if user.is_superuser:
            self.collection.update_one(
                {"username": username},
                {"$set": {
                    "is_superuser": False
                }}
            )
        else:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="remove superuser only permitted by superusers",
                headers={"WWW-Authenticate": "Bearer"},
            )

        return Response(content=f"Successfully reset the password for {username}!")

    def delete(
        self,
        username: str,
        token: Annotated[str, Depends(oauth2_scheme)]
    ) -> Response:
        user = self.get_current_user(token)
        if user.username == username or user.is_superuser:
            self.collection.delete_one({"username": username})
        else:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="delete operation only permitted by the owner or superusers",
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
