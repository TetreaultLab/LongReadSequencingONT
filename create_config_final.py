import argparse
import toml

parser = argparse.ArgumentParser(
    prog="PipelineLongReads",
    description="LongReadSequencing pipeline.",
    )

parser.add_argument(
    "config", type=str, help="Project config file, including path."
    )

args = parser.parse_args()

with open(args.config, "r") as f:
    toml_config = toml.load(f)
    print(toml_config)
