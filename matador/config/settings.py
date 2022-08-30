from collections import defaultdict


class Settings:
    """Light wrapper for global matador settings."""

    def __init__(self, settings=None):
        self.settings = defaultdict(dict)
        self.set = False

        if settings is not None:
            for key in settings:
                self.settings[key] = settings[key]
            self.set = True

    def __getitem__(self, key):
        # print("__getitem__", key)
        return self.settings[key]

    def get(self, key, default=None):
        return self.settings.get(key, default)

    def __setitem__(self, key, val):
        if isinstance(val, dict):
            self.settings[key].update(val)
        else:
            self.settings[key] = val

    def __iter__(self):
        return iter(self.settings)

    def __repr__(self):
        return str(self.settings)

    def reset(self):
        self.settings = defaultdict(dict)
        self.set = False
