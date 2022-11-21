import json

from configs.base_simulation_config import BaseSimulationConfig


class TeslaValveSimulationConfig(BaseSimulationConfig):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._current_direction = kwargs.get("_current_direction", 'easy')

    @property
    def current_direction(self):
        return self._current_direction

    @current_direction.setter
    def current_direction(self, direction):
        self._current_direction = direction

    def to_json(self):
        return json.dumps(self.__dict__)

    def __str__(self):
        return json.dumps(self.__dict__)
