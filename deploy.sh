#!/usr/bin/env bash

################################################################################
#
#   ASKCOS Deployment Utilities
#
#    ~ To streamline deployment commands ~
#
################################################################################

set -e  # exit with nonzero exit code if anything fails

# Create environment variable files from examples if they don't exist
if [ ! -f ".env" ]; then
  cp .env.example .env
fi

# Get docker compose variables from .env
source .env

# Default argument values
BUYABLES=""
CHEMICALS=""
REACTIONS=""
RETRO_TEMPLATES=""
FORWARD_TEMPLATES=""
SELEC_REFERENCES=""
DB_DROP=""
DROP_INDEXES=false
LOCAL=false
BACKUP_DIR=""
IGNORE_DIFF=false

COMMANDS=""
while (( "$#" )); do
  case "$1" in
    -h|--help|help)
      usage
      exit
      ;;
    -f|--compose-file)
      if [ -z "$COMPOSE_FILE" ]; then
        COMPOSE_FILE=$2
      else
        COMPOSE_FILE=$COMPOSE_FILE:$2
      fi
      shift 2
      ;;
    -p|--project-name)
      COMPOSE_PROJECT_NAME=$2
      shift 2
      ;;
    -l|--local)
      LOCAL=true
      shift 1
      ;;
    -v|--version)
      VERSION_NUMBER=$2
      shift 2
      ;;
    -b|--buyables)
      BUYABLES=$2
      shift 2
      ;;
    -c|--chemicals)
      CHEMICALS=$2
      shift 2
      ;;
    -x|--reactions)
      REACTIONS=$2
      shift 2
      ;;
    -r|--retro-templates)
      RETRO_TEMPLATES=$2
      shift 2
      ;;
    -t|--forward-templates)
      FORWARD_TEMPLATES=$2
      shift 2
      ;;
    -e|--selec-references)
      SELEC_REFERENCES=$2
      shift 2
      ;;
    -a|--append)
      echo 'The -a,--append argument is no longer necessary. Documents are appended by default.'
      shift 1
      ;;
    -z|--drop)
      DB_DROP="--drop"
      shift 1
      ;;
    -i|--drop-indexes)
      DROP_INDEXES=true
      shift 1
      ;;
    -n|--ignore-diff)
      IGNORE_DIFF=true
      shift 1
      ;;
    -d|--backup-directory)
      BACKUP_DIR=$2
      shift 2
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*) # any other flag
      echo "Error: Unsupported flag $1" >&2  # print to stderr
      exit 1
      ;;
    *) # preserve positional arguments
      COMMANDS="$COMMANDS $1"
      shift
      ;;
  esac
done

# Set positional arguments in their proper place
eval set -- "$COMMANDS"

# Export variables needed by docker compose
export VERSION_NUMBER
export COMPOSE_FILE
export COMPOSE_PROJECT_NAME

# Define various functions
diff-env() {
  # V2 update: dropped customization from diff since Django is no longer used
  if [ "$IGNORE_DIFF" = "true" ]; then
    return 0
  fi

  output1="$(diff -u .env .env.example)" || true
  if [ -n "$output1" ]; then
    echo -e "\033[91m*** WARNING ***\033[00m"
    echo "Local .env file is different from .env.example"
    echo "This could be due to changes to the example or local changes."
    echo "Please review the diff to determine if any changes are necessary:"
    echo
    echo "$output1"
    echo
  fi

  if [ -n "$output1" ]; then
    echo "Local configuration files differ from examples (see above)! What would you like to do?"
    echo "  (c)ontinue without changes  (use -n flag to skip prompt and continue in the future)"
    echo "  (o)verwrite your local files with the examples  (the above diff(s) will be applied)"
    echo "  (s)top and make changes manually"
    read -rp '>>> ' response
    case "$response" in
      [Cc])
        echo "Continuing without changes."
        ;;
      [Oo])
        echo "Overwriting local files."
        cp .env.example .env
        source .env
        ;;
      [Ss])
        echo "Stopping."
        exit 1
        ;;
      *)
        echo "Unrecognized option. Stopping."
        exit 1
        ;;
    esac
  fi
}

copy-https-conf() {
  echo "Using https nginx configuration."
  cp nginx.https.conf nginx.conf
  echo
  # Create SSL
  if [ ! -f "askcos.ssl.cert" ]; then
    echo "Creating SSL certificates."
    openssl req -new -newkey rsa:4096 -days 3650 -nodes -x509 -subj "/C=US/ST=MA/L=BOS/O=askcos/CN=askcos.$RANDOM.com" -keyout askcos.ssl.key -out askcos.ssl.cert
    echo
  fi
}

set-db-defaults() {
  # Set default values for seeding database if values are not already defined
  BUYABLES=${BUYABLES:-default}
  RETRO_TEMPLATES=${RETRO_TEMPLATES:-default}
  FORWARD_TEMPLATES=${FORWARD_TEMPLATES:-default}
  CHEMICALS=${CHEMICALS:-default}
  REACTIONS=${REACTIONS:-default}
  SELEC_REFERENCES=${SELEC_REFERENCES:-default}
}

run-mongo-js() {
  # arg 1 is js command
  docker compose -f compose.yaml exec -T mongo \
    bash -c 'mongosh --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin ${MONGO_HOST}/askcos --quiet --eval '"'$1'"
}

mongoimport() {
  # arg 1 is collection name
  # arg 2 is file path
  # arg 3 is a flag to pass to docker compose exec, e.g. -d to detach
  docker compose -f compose.yaml exec -T $3 mongo \
    bash -c 'gunzip -c '$2' | mongoimport --host ${MONGO_HOST} --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin --db askcos --collection '$1' --type json --jsonArray --numInsertionWorkers 8'${DB_DROP}
}

mongoexport() {
  # arg 1 is collection name
  # arg 2 is file path
  # arg 3 is a flag to pass to docker compose exec, e.g. -d to detach
  docker compose -f compose.yaml exec -T $3 mongo \
    bash -c 'mongoexport --host ${MONGO_HOST} --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin --db askcos --collection '$1' --type json --jsonArray | gzip > '$2
}

index-db() {
  if [ "$DROP_INDEXES" = "true" ]; then
    echo "Dropping existing indexes in mongo database..."
    run-mongo-js 'db.buyables.dropIndexes()'
    run-mongo-js 'db.chemicals.dropIndexes()'
    run-mongo-js 'db.reactions.dropIndexes()'
    run-mongo-js 'db.retro_templates.dropIndexes()'
    run-mongo-js 'db.sites_refs.dropIndexes()'
  fi
  echo "Adding indexes to mongo database..."
  run-mongo-js 'db.buyables.createIndex({smiles: 1, source: 1})'
  run-mongo-js 'db.chemicals.createIndex({smiles: 1, template_set: 1})'
  run-mongo-js 'db.reactions.createIndex({reaction_id: 1, template_set: 1})'
  run-mongo-js 'db.retro_templates.createIndex({index: 1, template_set: 1})'
  run-mongo-js 'db.sites_refs.createIndexes([{index: 1}, {reactant: 1}])'
  echo "Indexing complete."
  echo
}

count-mongo-docs() {
  echo "Buyables collection:          $(run-mongo-js "db.buyables.estimatedDocumentCount({})" | tr -d '\r') / 280469 expected (default)"
  echo "Chemicals collection:         $(run-mongo-js "db.chemicals.estimatedDocumentCount({})" | tr -d '\r') / 19188359 expected (default)"
  echo "Reactions collection:         $(run-mongo-js "db.reactions.estimatedDocumentCount({})" | tr -d '\r') / 5263537 expected (default)"
  echo "Retro template collection:    $(run-mongo-js "db.retro_templates.estimatedDocumentCount({})" | tr -d '\r') / 612013 expected (default)"
  echo "Forward template collection:  $(run-mongo-js "db.forward_templates.estimatedDocumentCount({})" | tr -d '\r') / 17089 expected (default)"
  echo "Sites ref collection:         $(run-mongo-js "db.sites_refs.estimatedDocumentCount({})" | tr -d '\r') / 123 expected (default)"
}

seed-db() {
  if [ -z "$BUYABLES" ] && [ -z "$CHEMICALS" ] && [ -z "$REACTIONS" ] && [ -z "$RETRO_TEMPLATES" ] && [ -z "$FORWARD_TEMPLATES" ] && [ -z "$REFERENCES" ]; then
    echo "Nothing to seed!"
    echo "Example usages:"
    echo "    bash deploy.sh seed-db -r default                  seed only the default retro templates"
    echo "    bash deploy.sh seed-db -r templates.json.gz        seed retro templates from local file templates.json.gz"
    echo "    bash deploy.sh set-db-defaults seed-db             seed all default collections"
    return
  fi

  echo "Starting the mongo container for seeding db"
  docker compose -f compose.yaml up -d mongo
  sleep 3

  echo "Seeding mongo database..."
  DATA_DIR="/usr/local/askcos-data/db"

  if [ "$BUYABLES" = "default" ]; then
    echo "Loading default buyables data..."
    buyables_file="${DATA_DIR}/buyables/buyables.json.gz"
    mongoimport buyables "$buyables_file"
  elif [ -f "$BUYABLES" ]; then
    echo "Loading buyables data from $BUYABLES in background..."
    buyables_file="${DATA_DIR}/buyables/$(basename "$BUYABLES")"
    docker cp "$BUYABLES" "${COMPOSE_PROJECT_NAME}"_mongo_1:"$buyables_file"
    mongoimport buyables "$buyables_file"
  fi

  if [ "$CHEMICALS" = "default" ]; then
    echo "Loading default chemicals data..."
    chemicals_file="${DATA_DIR}/historian/chemicals.json.gz"
    mongoimport chemicals "$chemicals_file"
    chemicals_file="${DATA_DIR}/historian/historian.pistachio.json.gz"
    DB_DROP="" mongoimport chemicals "$chemicals_file"
    chemicals_file="${DATA_DIR}/historian/historian.bkms_metabolic.json.gz"
    DB_DROP="" mongoimport chemicals "$chemicals_file"
  elif [ "$CHEMICALS" = "pistachio" ]; then
    echo "Loading pistachio chemicals data in background..."
    chemicals_file="${DATA_DIR}/historian/historian.pistachio.json.gz"
    mongoimport chemicals "$chemicals_file"
  elif [ "$CHEMICALS" = "bkms" ]; then
    echo "Loading bkms chemicals data in background..."
    chemicals_file="${DATA_DIR}/historian/historian.bkms_metabolic.json.gz"
    mongoimport chemicals "$chemicals_file"
  elif [ -f "$CHEMICALS" ]; then
    echo "Loading chemicals data from $CHEMICALS in background..."
    chemicals_file="${DATA_DIR}/historian/$(basename "$CHEMICALS")"
    docker cp "$CHEMICALS" "${COMPOSE_PROJECT_NAME}"_mongo_1:"$chemicals_file"
    mongoimport chemicals "$chemicals_file"
  fi

  if [ "$REACTIONS" = "default" ]; then
    echo "Loading default reactions data..."
    reactions_file="${DATA_DIR}/historian/reactions.pistachio.json.gz"
    mongoimport reactions "$reactions_file"
    reactions_file="${DATA_DIR}/historian/reactions.cas.min.json.gz"
    DB_DROP="" mongoimport reactions "$reactions_file"
    reactions_file="${DATA_DIR}/historian/reactions.bkms_metabolic.json.gz"
    DB_DROP="" mongoimport reactions "$reactions_file"
  elif [ "$REACTIONS" = "pistachio" ]; then
    echo "Loading pistachio reactions data in background..."
    reactions_file="${DATA_DIR}/historian/reactions.pistachio.json.gz"
    mongoimport reactions "$reactions_file"
  elif [ "$REACTIONS" = "cas" ]; then
    echo "Loading cas reactions data in background..."
    reactions_file="${DATA_DIR}/historian/reactions.cas.min.json.gz"
    mongoimport reactions "$reactions_file"
  elif [ "$REACTIONS" = "bkms" ]; then
    echo "Loading bkms reactions data in background..."
    reactions_file="${DATA_DIR}/historian/reactions.bkms_metabolic.json.gz"
    mongoimport reactions "$reactions_file"
  elif [ -f "$REACTIONS" ]; then
    echo "Loading reactions data from $REACTIONS in background..."
    reactions_file="${DATA_DIR}/historian/$(basename "$REACTIONS")"
    docker cp "$REACTIONS" "${COMPOSE_PROJECT_NAME}"_mongo_1:"$reactions_file"
    mongoimport reactions "$reactions_file"
  fi

  if [ "$RETRO_TEMPLATES" = "default" ]; then
    echo "Loading default retrosynthetic templates..."
    retro_file="${DATA_DIR}/templates/retro.templates.reaxys.json.gz"
    mongoimport retro_templates "$retro_file"
    retro_file="${DATA_DIR}/templates/retro.templates.pistachio.json.gz"
    DB_DROP="" mongoimport retro_templates "$retro_file"
    retro_file="${DATA_DIR}/templates/retro.templates.cas.json.gz"
    DB_DROP="" mongoimport retro_templates "$retro_file"
    retro_file="${DATA_DIR}/templates/retro.templates.bkms_metabolic.json.gz"
    DB_DROP="" mongoimport retro_templates "$retro_file"
    retro_file="${DATA_DIR}/templates/retro.templates.pistachio_ringbreaker.json.gz"
    DB_DROP="" mongoimport retro_templates "$retro_file"
    retro_file="${DATA_DIR}/templates/retro.templates.reaxys_biocatalysis.json.gz"
    DB_DROP="" mongoimport retro_templates "$retro_file"
  elif [ "$RETRO_TEMPLATES" = "pistachio" ]; then
    echo "Loading pistachio retrosynthetic templates..."
    retro_file="${DATA_DIR}/templates/retro.templates.pistachio.json.gz"
    mongoimport retro_templates "$retro_file"
  elif [ "$RETRO_TEMPLATES" = "cas" ]; then
    echo "Loading cas retrosynthetic templates..."
    retro_file="${DATA_DIR}/templates/retro.templates.cas.json.gz"
    mongoimport retro_templates "$retro_file"
  elif [ "$RETRO_TEMPLATES" = "bkms" ]; then
    echo "Loading bkms retrosynthetic templates..."
    retro_file="${DATA_DIR}/templates/retro.templates.bkms_metabolic.json.gz"
    mongoimport retro_templates "$retro_file"
  elif [ "$RETRO_TEMPLATES" = "ringbreaker" ]; then
    echo "Loading ringbreaker retrosynthetic templates..."
    retro_file="${DATA_DIR}/templates/retro.templates.pistachio_ringbreaker.json.gz"
    mongoimport retro_templates "$retro_file"
  elif [ -f "$RETRO_TEMPLATES" ]; then
    echo "Loading retrosynthetic templates from $RETRO_TEMPLATES ..."
    retro_file="${DATA_DIR}/templates/$(basename "$RETRO_TEMPLATES")"
    docker cp "$RETRO_TEMPLATES" "${COMPOSE_PROJECT_NAME}"_mongo_1:"$retro_file"
    mongoimport retro_templates "$retro_file"
  fi

  if [ "$FORWARD_TEMPLATES" = "default" ]; then
    echo "Loading default forward templates..."
    forward_file="${DATA_DIR}/templates/forward.templates.json.gz"
    mongoimport forward_templates "$forward_file"
  elif [ -f "$FORWARD_TEMPLATES" ]; then
    echo "Loading forward templates from $FORWARD_TEMPLATES ..."
    forward_file="${DATA_DIR}/templates/$(basename "$FORWARD_TEMPLATES")"
    docker cp "$FORWARD_TEMPLATES" "${COMPOSE_PROJECT_NAME}"_mongo_1:"$forward_file"
    mongoimport forward_templates "$forward_file"
  fi

  if [ "$SELEC_REFERENCES" = "default" ]; then
    echo "Loading default model references..."
    refs_file="${DATA_DIR}/references/site_selectivity.refs.json.gz"
    mongoimport sites_refs "$refs_file"
  fi

  index-db
  echo "Seeding db completed."
  count-mongo-docs
  docker compose -f compose.yaml rm -sf mongo
  echo
}

generate-deployment-scripts-in-docker() {
  if [ -z "${ASKCOS_REGISTRY}" ]; then
    export ASKCOS_REGISTRY=registry.gitlab.com/mlpds_mit/askcosv2
  fi

  echo "Building image for askcos2_core, runtime: docker"
  docker build -f Dockerfile -t ${ASKCOS_REGISTRY}/askcos2_core:2.0 .

  docker run --rm \
    -e ASKCOS2_CORE_DIR="$PWD" \
    -v "${PWD%/*}":/ASKCOSv2 \
    -t ${ASKCOS_REGISTRY}/askcos2_core:2.0 \
    python scripts/pre_deploy.py
}

generate-deployment-scripts() {
  python scripts/pre_deploy.py
}

get-images() {
  echo "Getting images using the script from the latest deployment directory"
  echo "    bash deployment/deployment_latest/get_images.sh"
  bash deployment/deployment_latest/get_images.sh
  echo "Images ready."
}

download-db-data() {
  echo "Downloading data for db using scripts/download_data.sh"
  bash scripts/download_data.sh
  echo "Data downloaded."
}

start-services() {
  echo "Start services using the script from the latest deployment directory"
  echo "    bash deployment/deployment_latest/start_services.sh"
  bash deployment/deployment_latest/start_services.sh
  echo "Services started."
}

stop-services() {
  echo "Stop services using the script from the latest deployment directory"
  echo "    bash deployment/deployment_latest/stop_services.sh"
  bash deployment/deployment_latest/stop_services.sh
  echo "Services stopped."
}

remove-volumes() {
  docker compose -f compose.yaml down -v
}

export_volume() {
  volume=$1
  directory=$2
  filename=$3
  volume_name=${COMPOSE_PROJECT_NAME}_${volume}
  docker run --rm -v "${volume_name}":/src -v "${directory}":/dest alpine tar -czf /dest/"${filename}" /src
}

import_volume() {
  volume=$1
  directory=$2
  filename=$3
  volume_name=${COMPOSE_PROJECT_NAME}_${volume}
  docker run --rm -v "${volume_name}":/dest -v "${directory}":/src alpine tar -xzf /src/"${filename}" -C /dest --strip 1
}

backup() {
  if [ -z "$BACKUP_DIR" ]; then
    BACKUP_DIR="$(pwd)/backup/$(date +%Y%m%d%s)"
  fi
  mkdir -p "${BACKUP_DIR}"
  echo "Backing up data to ${BACKUP_DIR}"
  echo "This may take a few minutes..."
  export_volume mongo_data "${BACKUP_DIR}" mongo_data.tar.gz
  echo "Backup complete."
}

restore() {
  if [ -z "$BACKUP_DIR" ]; then
    BACKUP_DIR="$(pwd)/backup/$(ls -t backup | head -1)"
  fi
  echo "Restoring data from ${BACKUP_DIR}"
  echo "This may take a few minutes..."
  import_volume mongo_data "${BACKUP_DIR}" mongo_data.tar.gz
  echo "Restore complete."
}

post-update-message() {
  echo
  echo -e "\033[92m================================================================================\033[00m"
  echo
  echo "The local ASKCOS deployment has been updated to version ${VERSION_NUMBER}!"
  echo
  echo "Please note the following items which may require further action:"
  echo
  echo "1) ASKCOS 2022.01 added reference data for the CAS template relevance model."
  echo "   If you have not done so already, you should import the required data:"
  echo
  echo "       bash deploy.sh seed-db -x cas"
  echo
  echo "2) ASKCOS 2021.10 added two new template relevance models: one for enzymatic"
  echo "   reactions from BKMS and one for ring-breaking reactions from Pistachio."
  echo "   If you have not done so already, you should import the required data:"
  echo
  echo "       bash deploy.sh seed-db -r bkms -c bkms -x bkms"
  echo "       bash deploy.sh seed-db -r ringbreaker"
  echo
  echo "3) ASKCOS 2021.07 added a new SciFinder/CAS-based template relevance model."
  echo "   If you have not done so already, you should import the required data:"
  echo
  echo "       bash deploy.sh seed-db -r cas"
  echo
  echo "   Please note that chemical historian data is not provided at this time."
  echo
  echo "4) ASKCOS 2021.04 added new reference data for the site selectivity model"
  echo "   and the Pistachio template relevance model."
  echo "   If you have not done so already, you should import the data as follows:"
  echo
  echo "       bash deploy.sh seed-db -e default -x default"
  echo
  echo "5) ASKCOS 2023.01 added a new template relevance model trained on"
  echo "   a Reaxys Enzymatic template set."
  echo "   If you have not done so already, you should import the data as follows:"
  echo "       bash deploy.sh seed-db -r reaxys_enzymatic"
  echo
  echo "For reference, past messages can be viewed using 'bash deploy.sh old-messages'."
  echo
  echo "                      ~~~ Thank you for using ASKCOS! ~~~"
  echo
  echo -e "\033[92m================================================================================\033[00m"
}

old-messages() {
  echo
  echo -e "\033[92m================================================================================\033[00m"
  echo
  echo "Post-update notes from past releases:"
  echo
  echo "1) ASKCOS 2020.07 changed the default MongoDB index types for much faster look-ups."
  echo "   If you have not done so already, you should recreate the MongoDB indexes:"
  echo
  echo "       bash deploy.sh index-db --drop-indexes"
  echo
  echo "2) ASKCOS 2020.10 added updated buyables data with more sources."
  echo "   If you have not done so already, you can import the data as follows:"
  echo
  echo "       bash deploy.sh seed-db -b default"
  echo
  echo "   If no custom buyables have been added and you would like to remove existing"
  echo "   buyables data before importing, you can add the '--drop' argument."
  echo
  echo "3) ASKCOS 2020.10 added a new Pistachio-based template relevance model."
  echo "   If you have not done so already, you should import the required data:"
  echo
  echo "       bash deploy.sh seed-db -c pistachio -r pistachio"
  echo
  echo "4) ASKCOS 2021.04 added new reference data for the site selectivity model"
  echo "   and the Pistachio template relevance model."
  echo "   If you have not done so already, you should import the data as follows:"
  echo
  echo -e "\033[92m================================================================================\033[00m"

}

# Handle positional arguments, which should be commands
if [ $# -eq 0 ]; then
  # No arguments
  echo "Must provide a valid task, e.g. deploy|update."
  echo "See 'deploy.sh help' for more options."
  exit 1;
else
  for arg in "$@"
  do
    case "$arg" in
      clean-data | start-db-services | save-db | seed-db | copy-https-conf | pull-images | \
      start-web-services | start-ml-servers | start-celery-workers | set-db-defaults | count-mongo-docs | \
      backup | restore | index-db | diff-env | post-update-message | old-messages | start-cpp-treebuilder-experimental )
        # This is a defined function, so execute it
        $arg
        ;;
      pre-deploy)
        copy-https-conf
        diff-env
        generate-deployment-scripts
        get-images
        download-db-data
        set-db-defaults
        seed-db
        ;;
      deploy)
        # Normal first deployment, do everything (pre-deploy + start-backend-services)
        copy-https-conf
        diff-env
        generate-deployment-scripts
        get-images
        download-db-data
        set-db-defaults
        seed-db
        start-services
        ;;
      update)
        # Update an existing configuration, database seeding is not performed
        copy-https-conf
        diff-env
        generate-deployment-scripts
        get-images
        start-services
        post-update-message
        ;;
      start)
        # (Re)start existing deployment
        start-services
        ;;
      stop)
        # Stop and remove currently running containers
        stop-services
        ;;
      clean)
        # Clean up current deployment including all data volumes
        echo "This will stop and remove all containers and also remove all data volumes, including user data."
        echo "Are you sure you want to continue?"
        read -rp 'Continue (y/N): ' response
        case "$response" in
          [Yy] | [Yy][Ee][Ss])
            echo "Cleaning deployment."
            stop-services
            remove-volumes
            ;;
          *)
            echo "Doing nothing."
            ;;
        esac
        ;;
      restart)
        # Stop and remove currently running containers before starting
        stop-services
        start-services
        ;;
      *)
        echo "Error: Unsupported command $arg" >&2  # print to stderr
        exit 1;
    esac
  done
fi
