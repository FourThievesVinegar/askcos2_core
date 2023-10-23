from datetime import datetime, timedelta
from fastapi import Depends, HTTPException, status
from fastapi.security import OAuth2PasswordBearer, OAuth2PasswordRequestForm
from jose import jwt
from passlib.context import CryptContext
from pydantic import BaseModel
from typing import Annotated
from utils.registry import get_util_registry

SECRET_KEY = "1d21d866fe6c71df2e2ab948e3d7d6bcf369c3d0f1f99d71570dc48365856c49"
ALGORITHM = "HS256"
# ACCESS_TOKEN_EXPIRE_MINUTES = 60
ACCESS_TOKEN_EXPIRE_MINUTES = 43200

# The scheme itself returns a str (i.e., the token) when called
# OAuth2PasswordBearer seems to have handled the conversion of token_dict -> token
oauth2_scheme = OAuth2PasswordBearer(tokenUrl="api/admin/token")
# optional_oauth2_scheme = OAuth2PasswordBearer(
#     tokenUrl="api/admin/token",
#     auto_error=False
# )

pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")

credentials_exception = HTTPException(
    status_code=status.HTTP_401_UNAUTHORIZED,
    detail="Could not validate credentials",
    headers={"WWW-Authenticate": "Bearer"},
)


class Token(BaseModel):
    access_token: str
    token_type: str = "bearer"


def create_access_token(data: dict, expires_delta: timedelta | None = None) -> str:
    to_encode = data.copy()
    if expires_delta:
        expire = datetime.utcnow() + expires_delta
    else:
        expire = datetime.utcnow() + timedelta(minutes=15)
    to_encode.update({"exp": expire})
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)

    return encoded_jwt


# The login function is directly bound to the token endpoint in app.py
# Note that the response_model is Token (auto-converted from token_dict)
async def login_for_access_token(
    form_data: Annotated[OAuth2PasswordRequestForm, Depends()]
) -> dict[str, str]:
    user_controller = get_util_registry().get_util(module="user_controller")
    user = user_controller.get_user_by_name(username=form_data.username)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username",
            headers={"WWW-Authenticate": "Bearer"},
        )
    if not pwd_context.verify(form_data.password, user.hashed_password):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect password",
            headers={"WWW-Authenticate": "Bearer"},
        )
    access_token_expires = timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
    access_token = create_access_token(
        data={"sub": user.username},
        expires_delta=access_token_expires
    )
    token_dict = {
        "access_token": access_token,
        "token_type": "bearer"
    }

    return token_dict
