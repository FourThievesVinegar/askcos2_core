import argparse
import importlib
import inspect
import os
import shutil
import subprocess
import sys
import yaml
from datetime import datetime

sys.path.append("")         # critical for importlib to work


def get_parser():
    parser = argparse.ArgumentParser("clone_repos")
    parser.add_argument("--config", help="Path to the config file",
                        type=str, default="configs/module_config_full.py")

    return parser


def run_and_printchar(args, cwd="."):
    # magic to pass subprocess printout to stdout. DO NOT CHANGE
    try:
        pipe = subprocess.Popen(
            args,
            bufsize=0,
            shell=False,
            cwd=cwd,
            stdout=None,
            stderr=subprocess.PIPE)
    except Exception as e:
        print(f"{inspect.stack()[1][3]} error: {str(e)}")
        return False

    started = False
    for line in pipe.stderr:
        line = line.decode("utf-8", "replace")
        if started:
            splited = line.split()
            if len(splited) == 9:
                percentage = splited[6]
                speed = splited[7]
                remaining = splited[8]
                print(f"Downloaded {percentage} with {speed} per second "
                      f"and {remaining} left.", flush=True, end="\r")
        elif line == os.linesep:
            started = True


def main():
    assert os.path.basename(os.path.abspath("..")) == "ASKCOSv2", \
        "Please ensure the askcos2_core is under the directory ASKCOSv2."

    if not os.path.exists(".env"):
        shutil.copy2(src=".env.example", dst=".env")

    args, _ = get_parser().parse_known_args()
    if args.config.endswith(".py"):
        m = args.config[:-3]
        m = m.replace("/", ".")
        m = importlib.import_module(m)
        module_config = getattr(m, "module_config")
    else:
        raise ValueError(f"Only .py configs are supported!")
    # module_config = validate(module_config) # TODO

    # set up the deployment directory
    cwd = os.environ.get("ASKCOS2_CORE_DIR", os.getcwd())
    os.makedirs("./deployment", exist_ok=True)
    dt = datetime.strftime(datetime.now(), "%y%m%d-%H%Mh")
    deployment_dir = f"./deployment/deployment_{dt}"

    os.makedirs(deployment_dir, exist_ok=True)
    if os.path.exists("./deployment/deployment_latest"):
        os.remove("./deployment/deployment_latest")
    os.symlink(os.path.basename(deployment_dir), "./deployment/deployment_latest")
    shutil.copy2(args.config, deployment_dir)

    runtime = module_config["global"]["container_runtime"]
    enable_gpu = module_config["global"]["enable_gpu"]
    image_policy = module_config["global"]["image_policy"]

    cmds_get_images = ["#!/bin/bash\n", "\n", "source .env\n", "\n"]
    cmds_start_services = ["#!/bin/bash\n", "\n", "source .env\n", "\n"]
    cmds_stop_services = ["#!/bin/bash\n", "\n", "source .env\n", "\n"]

    if runtime == "docker":
        cmds = [
            "echo Starting db services\n",
            "docker compose -f compose.yaml up -d rabbitmq redis mongo\n",
            "sleep 5\n"
            "\n"
        ]
        cmds_start_services.extend(cmds)

    for module, to_start in module_config["modules_to_start"].items():
        if not to_start:
            continue

        # clone module-to-start, if not existing
        repo = module_config[module]["repo"]
        output_dir = repo.split("askcosv2/")[-1][:-4]
        output_dir = f"../{output_dir}"

        is_newly_cloned = False
        if os.path.exists(output_dir):
            print(f"{output_dir} already exists, skipping cloning.")
            command = ["git", "pull"]
            run_and_printchar(command, cwd=output_dir)
        else:
            os.makedirs(output_dir)
            print(f"Cloning {repo} into {output_dir}")
            command = ["git", "clone", repo, output_dir]
            run_and_printchar(command)
            is_newly_cloned = True

        # load deployment config per repo, if specified
        # FIXME: this assumes "deployment.yaml" does exist, which should be the case
        try:
            config_file = module_config[module]["deployment"]["deployment_config"]
        except KeyError:
            raise ValueError(f"Deployment_config not specified for {module}")

        config_file = os.path.join(output_dir, config_file)
        with open(config_file, "r") as f:
            deployment_config = yaml.safe_load(f)

        # download module data, if any
        if is_newly_cloned:
            download_command = deployment_config["commands"].get("download", None)
            if download_command:
                print(f"Download command found, downloading data for {module}")
                run_and_printchar(download_command.split(), cwd=output_dir)
            else:
                print(f"No download command found for {module}, skipping")

        # prepare commands for getting the images
        if enable_gpu and module_config[module]["deployment"].get("use_gpu"):
            device = "gpu"
        else:
            device = "cpu"

        if image_policy == "build_all":
            if module in ["legacy_descriptors", "legacy_solubility"]:
                # only two legacy black-box images remaining
                # Hardcoded, DO NOT CHANGE
                cmd = str(deployment_config[runtime][device]["image"])
                cmd = f"docker pull {cmd}"
                cmds = [
                    f"echo Pulling image for module: {module}, "
                    f"runtime: {runtime}, device: {device}\n",
                    f"cd {output_dir}\n",
                    f"{cmd}\n",
                    f"cd {cwd}\n",
                    "\n"
                ]
            else:
                cmd = str(deployment_config[runtime][device]["build"])
                cmds = [
                    f"echo Building image for module: {module}, "
                    f"runtime: {runtime}, device: {device}\n",
                    f"cd {output_dir}\n",
                    f"{cmd}\n",
                    f"cd {cwd}\n",
                    "\n"
                ]
            cmds_get_images.extend(cmds)
        else:
            raise NotImplementedError

        # prepare commands for starting and stopping services
        cmd = str(deployment_config[runtime][device]["start"])
        cmds = [
            f"echo Starting service for module: {module}, "
            f"runtime: {runtime}, device: {device}\n",
            f"cd {output_dir}\n",
            f"{cmd}\n",
            f"cd {cwd}\n",
            "\n"
        ]

        cmds_start_services.extend(cmds)

        cmd = str(deployment_config["commands"][f"stop-{runtime}"])
        cmds = [
            f"echo Stopping service for module: {module}, runtime: {runtime}\n",
            f"cd {output_dir}\n",
            f"{cmd}\n",
            f"cd {cwd}\n",
            "\n"
        ]
        cmds_stop_services.extend(cmds)

    # extend with commands for app and celery
    with open("deployment.yaml", "r") as f:
        core_deployment_config = yaml.safe_load(f)

    if image_policy == "build_all":
        cmd = str(core_deployment_config[runtime]["cpu"]["build"])
        cmds = [
            f"echo Building image for askcos2_core, runtime: {runtime}\n",
            f"{cmd}\n",
            "\n"
        ]
        cmds_get_images.extend(cmds)
    else:
        raise NotImplementedError

    cmd = str(core_deployment_config[runtime]["cpu"]["start"])
    cmds = [
        f"echo Starting services for askcos2_core app and celery_workers, "
        f"runtime: {runtime}\n",
        f"{cmd}\n",
        "\n"
    ]
    cmds_start_services.extend(cmds)

    cmd = str(core_deployment_config["commands"][f"stop-{runtime}"])
    cmds = [
        f"echo Stopping services for db, askcos2_core app and any celery workers, "
        f"runtime: {runtime}\n",
        f"{cmd}\n",
        "\n"
    ]
    cmds_stop_services.extend(cmds)

    print(f"Writing deployment commands into {deployment_dir}")
    ofn = os.path.join(deployment_dir, "get_images.sh")
    with open(ofn, "w") as of:
        of.writelines(cmds_get_images)

    ofn = os.path.join(deployment_dir, "start_services.sh")
    with open(ofn, "w") as of:
        of.writelines(cmds_start_services)

    ofn = os.path.join(deployment_dir, "stop_services.sh")
    with open(ofn, "w") as of:
        of.writelines(cmds_stop_services)


if __name__ == "__main__":
    main()
