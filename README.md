# askcos2-core

The core FastAPI app (a.k.a. the API gateway) for ASKCOSv2. The functional modules have their own repositories, which will be automatically pulled during setup based on the provided config file.

## Running ASKCOSv2 backend

In the general sense, the ASKCOSv2 backend consists of the API gateway and the microservices that serve the functional modules. We therefore need to set up the microservices for those modules, before we can call them via the API gateway.

### Option 1: the easy way
We have made use of Makefile and the deployment script to make deployment with the default config easy. First make sure that `askcos2_core` is under a folder named `ASKCOSv2`. Then while in `askcos2_core`,
```
make deploy
```
which will run the pre-deployment setup by cloning the backend repos and downloading the data. This command also builds all the docker images locally and seeds the mongo database with all provided data.

To temporarily stop a deployment
```
make stop
```

To (re)start a deployment
```
make start
```

### Option 2: manual deployment with more customization
#### Step 1. Configure and set up the module repositories

First make sure that `askcos2_core` is under a folder named `ASKCOSv2`. Then while in `askcos2_core`, pick a configuration file under `./configs` for your use case. Use `module_config_full.py` if you are not sure. Advanced users can create own config based on the pattern in `module_config_full.py`. Then prepare for deployment with the config file, for example, by
```
python scripts/pre_deploy.py --config configs/module_config_full.py
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

The pre-deployment script then generate the deployment scripts from templates given the config, and put them under a newly created directory `./deployment/deployment_TIMESTAMP/`, which is automatically aliased as `./deployment/deployment_latest/`. This helps keep track of any past deployment.

#### Step 2. Get docker/singularity images

```
bash deployment/deployment_latest/get_images.sh
```

#### Step 3. Start services

To start the ASKCOSv2 backend locally, run
```
bash deployment/deployment_latest/start_services.sh
```
which will start the backend microservices, the celery workers (not supported for singularity), and the core FastAPI app. Note that the services are started in the background (i.e., in detached mode). So they would need to be explicitly stopped if no longer in use, by
```
bash deployment/deployment_latest/stop_services.sh
```
For more advanced deployment with K8s and/or helm, please refer to the repo `askcos2_deploy`.

## Sample API usage

The backend microservices are managed centrally through the FastAPI app, running at http://0.0.0.0:9100 by default. To use any endpoint, you can first look at the API definition at http://0.0.0.0:9100/docs. You should mostly be making the synchronous call, unless you are interacting with the frontend. The POST requests with desired data format defined in docs should be sent to the `call_sync` endpoint, for example, to http://0.0.0.0:9100/site_selectivity/call_sync.

## Unit Test
With the services started, run
```
V2_HOST=0.0.0.0 make test
```