import copy
import unittest
import os
import shutil
from matador.config import load_custom_settings, set_settings
from matador.config.config import DEFAULT_SETTINGS

DUMMY_SETTINGS = {
    "mongo": {"host": "blah", "port": 666},
    "plotting": {"style": "matador"},
    "this_is_a_test": {True: "it is"},
    "run3": {
        "castep_executable": "castep",
        "optados_executable": "optados",
        "scratch_prefix": ".",
    },
}

REAL_PATH = "/".join(__file__.split("/")[:-1]) + "/"


class ConfigTest(unittest.TestCase):
    """Test config loading."""

    def tearDown(self):
        from matador.config import SETTINGS

        SETTINGS.reset()

    def setUp(self):
        from matador.config import SETTINGS

        SETTINGS.reset()

    def testLoadNamedCustomSettings(self):
        """Test custom config."""
        settings = load_custom_settings(
            config_fname=(REAL_PATH + "data/custom_config.yml"), no_quickstart=True
        )
        self.assertEqual(settings.settings, DUMMY_SETTINGS)
        from matador.config import SETTINGS

        self.assertEqual(SETTINGS.settings, DUMMY_SETTINGS)
        SETTINGS.reset()

    def testSetSettings(self):
        set_settings(DUMMY_SETTINGS)
        from matador.config import SETTINGS

        self.assertEqual(SETTINGS.settings, DUMMY_SETTINGS)

        SETTINGS["backend"] = "mongo"
        self.assertEqual(SETTINGS.settings["backend"], "mongo")
        SETTINGS.reset()

    def testLoadUserDefaultSettings(self):
        """Test default config."""
        exists = False
        try:
            if os.path.isfile(os.path.expanduser("~/.matadorrc")):
                exists = True
                shutil.copy(
                    os.path.expanduser("~/.matadorrc"),
                    os.path.expanduser("~/.matadorrc_bak"),
                )
            shutil.copy(
                REAL_PATH + "data/custom_config.yml", os.path.expanduser("~/.matadorrc")
            )
            settings = load_custom_settings(no_quickstart=True)
            self.assertEqual(settings.settings, DUMMY_SETTINGS)
            os.remove(os.path.expanduser("~/.matadorrc"))
            if exists:
                shutil.copy(
                    os.path.expanduser("~/.matadorrc_bak"),
                    os.path.expanduser("~/.matadorrc"),
                )
                os.remove(os.path.expanduser("~/.matadorrc_bak"))
        except Exception as oops:
            if exists:
                shutil.copy(
                    os.path.expanduser("~/.matadorrc_bak"),
                    os.path.expanduser("~/.matadorrc"),
                )
                os.remove(os.path.expanduser("~/.matadorrc_bak"))
            raise oops

    def testLoadDefaultSettings(self):
        """Test default config."""
        settings = load_custom_settings(
            config_fname="definitely_doesnt_exist.yml", no_quickstart=True
        )
        self.assertEqual(settings.settings, DEFAULT_SETTINGS)

        settings.reset()

    def testSetDefaultSettings(self):
        """Test default config."""
        set_settings(DEFAULT_SETTINGS)
        from matador.config import SETTINGS

        new_settings = copy.deepcopy(DUMMY_SETTINGS)
        new_settings["mongo"].pop("host")

        set_settings(new_settings)
        self.assertEqual(
            SETTINGS.settings["mongo"]["host"], DEFAULT_SETTINGS["mongo"]["host"]
        )

        SETTINGS.reset()
