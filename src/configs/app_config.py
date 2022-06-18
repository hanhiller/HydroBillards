from pathlib import Path


class AppConfig:
    """
    Application specific configuration. Here application means a part of this program that does not
    directly deal with simulation business logic.
    """
    def __init__(self):
        self._outputPath = Path.cwd()
