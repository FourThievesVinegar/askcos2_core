import uvicorn
from fastapi import APIRouter, FastAPI
from wrappers.registry import get_wrapper_registry

app = FastAPI()

for wrapper in get_wrapper_registry():
    for prefix in wrapper.prefixes:
        router = APIRouter(prefix=prefix)
        # Bind specified wrapper.method to urls
        # /{wrapper.prefixes}/method_name
        for method_name in wrapper.methods_to_bind:
            router.add_api_route(
                path=f"/{method_name}",
                endpoint=getattr(wrapper, method_name),
                methods=["POST"]
            )
        app.include_router(router)

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=9100)
