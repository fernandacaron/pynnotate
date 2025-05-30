import argparse

def load_config(config_path):
    
    import yaml

    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def merge_args_with_config(config: dict, args: argparse.Namespace, parser: argparse.ArgumentParser) -> dict:
    final = config.copy()

    for key in vars(args):
        arg_val = getattr(args, key)
        default_val = parser.get_default(key)

        if arg_val != default_val and arg_val is not None:
            final[key] = arg_val

    return final