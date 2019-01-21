from enum import Enum


class Components(Enum):
    OIL = 0
    WATER = 1


class FlowComponents(Enum):
    minus = 0
    plus = 1


class Dim(Enum):
    i = 0
    j = 1
    k = 2


class FlowCell(Enum):
    left_cell = 'left_cell'
    right_cell = 'right_cell'
