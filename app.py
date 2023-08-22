import uvicorn
from fastapi import APIRouter, FastAPI
from utils.registry import get_util_registry
from wrappers.registry import get_wrapper_registry

util_registry = get_util_registry()
wrapper_registry = get_wrapper_registry()

app = FastAPI()
router = APIRouter(prefix="/api/admin")
router.add_api_route(
    path="/get_backend_status",
    endpoint=wrapper_registry.get_backend_status,
    methods=["GET"]
)
app.include_router(router)

for wrapper in wrapper_registry:
    for prefix in wrapper.prefixes:
        router = APIRouter(prefix=f"/api/{prefix}")
        # Bind specified wrapper.method to urls /{wrapper.prefixes}/method_name
        for method_name, bind_types in wrapper.methods_to_bind.items():
            router.add_api_route(
                path=f"/{method_name}",
                endpoint=getattr(wrapper, method_name),
                methods=bind_types
            )
        app.include_router(router)

for util in util_registry:
    for prefix in util.prefixes:
        router = APIRouter(prefix=f"/api/{prefix}")
        # Bind specified util.method to urls /{util.prefixes}/method_name
        for method_name, bind_types in util.methods_to_bind.items():
            router.add_api_route(
                path=f"/{method_name}",
                endpoint=getattr(util, method_name),
                methods=bind_types
            )
        app.include_router(router)

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=9100)
