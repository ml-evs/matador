#!/usr/bin/env python
import unittest
import os
import shutil
from matador.config import load_custom_settings
from matador.config.config import DEFAULT_SETTINGS

DUMMY_SETTINGS = {'mongo': {'host': 'blah', 'port': 666}, 'plotting': {'style': 'matador'}, 'this_is_a_test': {True: 'it is'}}
REAL_PATH = '/'.join(__file__.split('/')[:-1]) + '/'


class ConfigTest(unittest.TestCase):
    """ Test config loading. """
    def testLoadNamedCustomSettings(self):
        """ Test custom config. """
        settings = load_custom_settings(config_fname=(REAL_PATH+'data/custom_config.yml'), override=True)
        self.assertEqual(settings, DUMMY_SETTINGS)

    def testLoadUserDefaultSettings(self):
        """ Test default config. """
        exists = False
        try:
            if os.path.isfile(os.path.expanduser('~/.matadorrc')):
                exists = True
                shutil.copy(os.path.expanduser('~/.matadorrc'), os.path.expanduser('~/.matadorrc_bak'))
            shutil.copy(REAL_PATH + 'data/custom_config.yml', os.path.expanduser('~/.matadorrc'))
            settings = load_custom_settings(override=True)
            self.assertEqual(settings, DUMMY_SETTINGS)
            os.remove(os.path.expanduser('~/.matadorrc'))
            if exists:
                shutil.copy(os.path.expanduser('~/.matadorrc_bak'), os.path.expanduser('~/.matadorrc'))
                os.remove(os.path.expanduser('~/.matadorrc_bak'))
        except Exception as oops:
            if exists:
                shutil.copy(os.path.expanduser('~/.matadorrc_bak'), os.path.expanduser('~/.matadorrc'))
                os.remove(os.path.expanduser('~/.matadorrc_bak'))
            raise oops

    def testLoadDefaultSettings(self):
        """ Test default config. """
        settings = load_custom_settings(config_fname='definitely_doesnt_exist.yml', override=True)
        self.assertEqual(settings, DEFAULT_SETTINGS)


if __name__ == '__main__':
    unittest.main()
