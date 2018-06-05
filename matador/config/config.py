# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the loading custom settings for plotting and
for database collections.

"""


import os
from traceback import print_exc
from matador.utils.print_utils import print_warning


DEFAULT_SETTINGS = {
    "mongo": {
        "db": "crystals",
        "host": "localhost",
        "port": 27017,
        "default_collection": "repo",
    },
    "plotting": {
        "default_style": "matador",
        "element_colours": "vesta_elements.ini",
    }
}


def load_custom_settings(config_fname=None, debug=False):
    """ Load mongodb settings dict from file given by fname, or from
    defaults. Hierarchy of filenames to check:

        1. ~/.matadorrc
        2. ~/.config/matadorrc

    Keyword Arguments:
        fname (str): filename of custom settings file.
        debug (bool): print settings on loading.

    """
    import yaml

    if config_fname is None:
        trial_user_fnames = [
            os.path.expanduser("~/.matadorrc"),
            os.path.expanduser("~/.config/matadorrc"),
            '/'.join(__file__.split('/')[:-1]) + 'matadorrc.yml'
        ]
        if sum([os.path.isfile(fname) for fname in trial_user_fnames]) > 1:
            print_warning("Found multiple user config files {}"
                          .format([val for val in trial_user_fnames if os.path.isfile(val)]))

        # use first file that exists in hierarchy
        for fname in trial_user_fnames:
            if os.path.isfile(fname):
                config_fname = fname
                break
        else:
            # otherwise load default
            config_fname = (
                "/".join(__file__.split("/")[:-1]) + "/matadorrc.yml"
            )

        print("Loading settings from {}".format(config_fname))

    custom_settings = {}
    if os.path.isfile(config_fname):
        try:
            with open(config_fname, 'r') as f:
                custom_settings = yaml.load(f)
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

    if debug:
        import json
        print(json.dumps(settings, indent=2))

    return settings
