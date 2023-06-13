import argparse
import importlib
import inspect
import os
import subprocess
import yaml


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
                      f"and {remaining} left.")
        elif line == os.linesep:
            started = True


def main():
    assert os.path.basename(os.path.abspath("..")) == "ASKCOSv2", \
        "Please ensure the askcos2_core is under the directory ASKCOSv2."

    args, _ = get_parser().parse_known_args()

    if args.config.endswith(".py"):
        m = args.config[:-3]
        m = m.replace("/", ".")
        m = importlib.import_module(m)
        module_config = getattr(m, "module_config")
    else:
        raise ValueError(f"Only .py configs are supported!")

    for module, to_start in module_config["modules_to_start"].items():
        if not to_start:
            continue

        # clone module-to-start, if not existing
        repo = module_config[module]["repo"]
        output_dir = repo.split("askcosv2/")[-1][:-4]
        output_dir = f"../{output_dir}"

        if os.path.exists(output_dir):
            print(f"{output_dir} already exists, skipping cloning.")
            continue

        os.makedirs(output_dir)
        print(f"Cloning {repo} into {output_dir}")
        command = ["git", "clone", repo, output_dir]
        run_and_printchar(command)

        # download module data, if any
        try:
            config_file = module_config[module]["deployment"]["deployment_config"]
        except KeyError:
            print(f"deployment_config not specified for {module}, skipping")
            continue

        config_file = os.path.join(output_dir, config_file)
        with open(config_file, "r") as f:
            deployment_config = yaml.safe_load(f)

        try:
            download_command = deployment_config["commands"]["download"]
        except KeyError:
            print(f"No download command found for {module}, skipping")
            continue

        print(f"Download command found, downloading data for {module}")
        run_and_printchar(download_command.split(), cwd=output_dir)


if __name__ == "__main__":
    main()
