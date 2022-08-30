# coding: utf-8
# Distributed under the terms of the MIT License.

""" The config simply loads the default or user-specified config. """


__all__ = ["load_custom_settings", "set_settings", "SETTINGS"]
__author__ = "Matthew Evans"
__maintainer__ = "Matthew Evans"

from .config import load_custom_settings, SETTINGS, set_settings
