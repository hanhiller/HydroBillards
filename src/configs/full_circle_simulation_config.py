from configs.base_simulation_config import BaseSimulationConfig
import json


class FullCircleSimulationConfig(BaseSimulationConfig):
    """
    Configuration to simulation full circle.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._constriction_width = kwargs.get("_constriction_width", 0.5)
        self._diameter = kwargs.get("_diameter", 30)
        self._inject_height = kwargs.get("_inject_height", 1.5)
        self._inject_width = kwargs.get("_inject_width", 0.6)

    @property
    def constriction_width(self):
        return self._constriction_width

    @constriction_width.setter
    def constriction_width(self, value):
        self._constriction_width = value

    @property
    def diameter(self):
        return self._diameter

    @diameter.setter
    def diameter(self, value):
        self._diameter = value

    @property
    def inject_height(self):
        return self._inject_height

    @inject_height.setter
    def inject_height(self, value):
        self._inject_height = value

    @property
    def inject_width(self):
        return self._inject_width

    @inject_width.setter
    def inject_width(self, value):
        self._inject_width = value

    def to_json(self):
        return json.dumps(self.__dict__)

    def __str__(self):
        return json.dumps(self.__dict__)
