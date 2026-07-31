'''
Created on 6 April 2026
read and write the YAML training recipes

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import json


#%% io
def load_config(path):
    path = Path(path)
    if path.suffix.lower() == '.json':
        with open(path, 'r') as f:
            return json.load(f)

    try:
        import yaml
    except ImportError as exc:
        raise ImportError('please install pyyaml to read yaml configs') from exc

    with open(path, 'r') as f:
        return yaml.safe_load(f)


def save_config(config, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if path.suffix.lower() == '.json':
        with open(path, 'w') as f:
            json.dump(config, f, indent=2)
        return

    try:
        import yaml
    except ImportError as exc:
        raise ImportError('please install pyyaml to write yaml configs') from exc

    with open(path, 'w') as f:
        yaml.safe_dump(config, f, sort_keys=False)


def get_section(config, section, defaults=None):
    values = dict(defaults or {})
    values.update(config.get(section, {}) or {})
    return values


def resolve_path(path, repo_root):
    path = Path(path)
    if path.is_absolute():
        return path
    return Path(repo_root) / path
