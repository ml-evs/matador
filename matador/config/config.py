# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the loading custom settings for plotting and
for database collections.

"""


import os
from pathlib import Path
from traceback import print_exc
from matador.utils.print_utils import print_warning
from matador.config.settings import Settings
from matador.config.quickstart import quickstart_settings

SETTINGS = Settings()

DEFAULT_SETTINGS = {
    "mongo": {
        "db": "crystals",
        "host": "localhost",
        "port": 27017,
        "default_collection": "repo",
    },
    "plotting": {
        "default_style": "matador",
        "element_colours": Path(__file__).parent.joinpath("vesta_elements.ini"),
    },
    "run3": {
        "scratch_prefix": ".",
        "castep_executable": "castep",
        "optados_executable": "optados",
    },
}

FILE_PATHS = {
    "mongo": ["default_collection_file_path", "scratch_prefix"],
    "plotting": ["element_colours"],
}


def load_custom_settings(
    config_fname=None, quiet=False, debug=False, no_quickstart=True
):
    """Load mongodb settings dict from file given by fname, or from
    defaults. Hierarchy of filenames to check:

        1. .matadorrc
        2. ~/.matadorrc
        3. ~/.config/matadorrc

    Keyword Arguments:
        fname (str): filename of custom settings file.
        quiet (bool): print nothing on loading.
        no_quickstart (bool): for testing purposes, override any questions.
        debug (bool): print settings on loading.

    """
    if SETTINGS.set:
        return SETTINGS

    import yaml

    if config_fname is None:
        trial_user_fnames = [
            Path(".matadorrc"),
            Path(os.path.expanduser("~/.matadorrc")),
            Path(os.path.expanduser("~/.config/matadorrc")),
            Path(__file__).parent.joinpath("matadorrc.yml"),
        ]

        # use first file that exists in hierarchy
        for fname in trial_user_fnames:
            if os.path.isfile(fname):
                config_fname = str(fname)
                break
        else:
            # otherwise offer to create one, or use default
            if not no_quickstart:
                config_fname = quickstart_settings()
            if config_fname is None:
                # otherwise load default
                print(
                    "No config provided: loading perhaps unsuitable defaults instead."
                )
                config_fname = (Path(__file__).parent.joinpath("matadorrc.yml"),)

        if not quiet:
            print("Loading settings from {}".format(config_fname))

    custom_settings = {}
    if os.path.isfile(config_fname):
        try:
            with open(config_fname, "r") as f:
                custom_settings = yaml.safe_load(f)
                if custom_settings is None:
                    custom_settings = {}
        except Exception:
            print_exc()
            raise SystemExit(
                "Failed to read custom settings file {} as YAML".format(config_fname)
            )
    else:
        print_warning(
            "Could not find {}, loading default settings...".format(config_fname)
        )

    settings = {}
    settings.update(DEFAULT_SETTINGS)
    settings.update(custom_settings)

    # check absolute/relative paths
    for key in FILE_PATHS:
        for setting in FILE_PATHS[key]:
            if settings.get(key) is not None:
                if settings[key].get(setting) is not None:
                    settings[key][setting] = str(
                        Path(os.path.expanduser(settings[key][setting])).resolve()
                    )

    if debug:
        import json

        print(json.dumps(settings, indent=2))

    set_settings(settings, override=False)

    return SETTINGS


def set_settings(settings, override=True):
    """Updates the global settings with the dictionary
    provided.

    Parameters:
        settings (dict): the key-value pairs to use.

    Keyword arguments:
        override (bool): whether or not to override
            previous settings, or just update.

    """

    if SETTINGS.set and not override:
        return

    for key in settings:
        if key not in SETTINGS:
            SETTINGS[key] = {}
        for _key in settings[key]:
            SETTINGS[key][_key] = settings[key][_key]

    SETTINGS.set = True
