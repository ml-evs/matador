from collections import defaultdict


class Settings:
    """ Light wrapper for global matador settings. """
    def __init__(self, settings=None):
        self.set = False
        if settings is None:
            self.settings = defaultdict(dict)
        else:
            self.settings = defaultdict(dict)
            for key in settings:
                self.settings[key] = settings[key]
            self.set = True

    def __getitem__(self, key):
        return self.settings[key]

    def __setitem__(self, key, val):
        self.settings[key] = val
