# askcos2-core
The microservice-based backend for ASKCOSv2. Note that the functional modules exist in their own repositories, which will be automatically pulled during setup. **askcos2-common** mainly houses the centralized FastAPI app, with some common utilities.

## Running ASKCOSv2 backend
### 1. Setup
#### 1.1 Configure and set up the module repositories

Pick a configuration file in **./configs** for your use case. Use **full.conf** if you are not sure. Advanced users can create own config based on the pattern in **full.conf**. Then clone the module repos with the config file, for example by
```
    sh scripts/clone_repos.sh --config configs/full.conf
```
You can check that required repos have been pulled under **./askcos2-core**.

#### 1.2. Pull the data files

We can now pull the big data files (databases, model checkpoints, etc.), again using the configuration file by
```
    sh scripts/pull_data.sh --config configs/full.conf
```

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

#### 1.4. Finalize
This step is for some final setup such as database migration.
```
    sh scripts/finalize.sh
```

### 2. Start services
To start services locally, run
```
    sh scripts/start_services.sh --config configs/full.conf
```
which relies on **docker-compose.yml**. The FastAPI app is able to figure out the dependencies automatically based on the config.

For more advanced deployment with K8s and/or helm, please refer to the deployment guide in **askcos2-deploy**.

### 3. Sample API usage
Once the services have been started, they are managed centrally through the FastAPI app, running at http://0.0.0.0:9100 by default. To use any endpoint, you can first look at the API definition at http://0.0.0.0:9100/docs. You should mostly be making the synchronous call, unless you are interacting with the front-end. The POST requests with desired data format defined in docs should be sent to the *call_sync* endpoint, for example, to http://0.0.0.0:9100/path_finder_mcts/call_sync.
