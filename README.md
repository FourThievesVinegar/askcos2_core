# askcos2-core
The core FastAPI app (a.k.a. the API gateway) for ASKCOSv2. The functional modules have their own repositories, which will be automatically pulled during setup based on the provided config file.

## Running ASKCOSv2 backend

In the general sense, the ASKCOSv2 backend consists of the API gateway and the microservices that serve the functional modules. We therefore need to set up the microservices for those modules, before we can call them via the API gateway.

### Step 1. Set up the microservices

#### 1.1 Configure and set up the module repositories

First make sure that `askcos2_core` is under a folder named `ASKCOSv2`. Then while in `askcos2_core`, pick a configuration file under `./configs` for your use case. Use `model_config_full.py` if you are not sure. Advanced users can create own config based on the pattern in `model_config_full.py`. Then clone the module repos with the config file, for example, by
```
    python scripts/clone_repos.py --config configs/model_config_full.py
```
The required repos specified in the config will also be put under `ASKCOSv2` automatically. A sample file organization will look like
```
ASKCOSv2
|-- askcos2_core
|-- site_selectivity
|-- forward_predictor
    |-- graph2smiles
    |-- augmented_transformer
...
```
If there are big data files used by the modules (as specified in the config file), they will also be downloaded and put into the correct places.

#### 1.3. Build/pull container images

Build all the docker images by
```
    sh scripts/build_all_dockers.sh --with_gpu_support
```
or if you need singularity images,
```
    sh scripts/build_all_singularities.sh --with_gpu_support
```
where the **--with_gpu_support** flag is optional. GPU-based images are recommended if you want to retrain the models, or if your use case needs high throughput.
You can also pull our prebuilt images with the corresponding scripts.

[//]: # (#### 1.4. Finalize)
[//]: # (This step is for some final setup such as database migration.)
[//]: # (```)
[//]: # (    sh scripts/finalize.sh)
[//]: # (```)

### 2. Start services
To start the ASKCOSv2 backend locally, run
```
    sh scripts/start_services.sh --config configs/model_config_full.py
```
which will start the backend microservices, the celery workers, and the core FastAPI app based on `compose.yaml` under this repo, as well as `deployment.yaml` under all module repos. For more advanced deployment with K8s and/or helm, please refer to the repo `askcos2_deploy`.

### 3. Sample API usage

The backend microservices are managed centrally through the FastAPI app, running at http://0.0.0.0:9100 by default. To use any endpoint, you can first look at the API definition at http://0.0.0.0:9100/docs. You should mostly be making the synchronous call, unless you are interacting with the frontend. The POST requests with desired data format defined in docs should be sent to the `call_sync` endpoint, for example, to http://0.0.0.0:9100/site_selectivity/call_sync.
