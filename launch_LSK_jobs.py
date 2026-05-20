import argparse
import os
import subprocess
import toml

parser = argparse.ArgumentParser(
    prog="Launch all LSK jobs in a LRS batch",
    description="change config and launch script",
)

parser.add_argument(
    "--email", type=str, required=True, help="email address for SLURM notification."
)
parser.add_argument("--analyses", type=str, required=True, help="which analyses to do.")
args = parser.parse_args()

email = args.email
analyses = args.analyses

cwd = os.getcwd()

subfolders = [f.path for f in os.scandir(cwd) if f.is_dir()]

# Loop over subfolders/samples
for s in subfolders:
    os.chdir(s)
    print("Current Working Directory:", os.getcwd())

    # Open config file and change email adress
    toml_config_initial = toml.load("config_initial.txt")
    toml_config_initial["general"]["email"] = email

    # Check if analyses are specified
    if toml_config_initial["general"]["analysis"] == [""]:
        toml_config_initial["general"]["analysis"] = analyses
    else:
        print(toml_config_initial["general"]["analysis"])

    # Save modification
    with open("config.toml", "w") as f:
        toml.dump(toml_config_initial, f)

    # Launch pipeline
    print("Launching pipeline")
    subprocess.run(
        [
            "python",
            "/lustre09/project/6019267/shared/tools/main_pipelines/long-read/LongReadSequencingONT/launch_pipeline.py",
            s + "/config.toml",
            "--test",
        ]
    )
