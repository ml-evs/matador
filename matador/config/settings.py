class Settings:
    def __init__(self, settings=None):
        self.set = False
        if settings is None:
            self.settings = {}
        else:
            self.settings = settings
            self.set = True

    def __getitem__(self, key):
        return self.settings[key]
