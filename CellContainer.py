import numpy as np

from .Layer import Layer
from .Cell import Cell


class CellContainer:
    def __init__(self):
        self.container = np.zeros((Layer.N_z, Layer.N_x, Layer.N_z), dtype=Cell)  # (z, x, y) x = rows, y - columns
        eq_index = 0
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    # TODO: Проверить действительно ли надо использовать такой номер уравнения
                    # Проверяем является ли клетка граничной, если да, то ставим, что она граничная
                    if (k == 0 or k == Layer.N_z) or(i == 0 or i == Layer.N_x) or (j == 0 or j == Layer.N_y):
                        self.container[k, i, j] = Cell(eq_index, True)
                    else:
                        self.container[k, i, j] = Cell(eq_index)
                    eq_index += 1

    def get_cells(self):
        return self.container

    def initialize_cell(self, cell):
        # TODO: initialize one cell
        for state in cell.get_cell_states():
            state.set_s_water(Layer.s_water_init * np.ones((Layer.dim, 1)))
            state.set_s_oil((1.0 - Layer.s_water_init) * np.ones((Layer.dim, 1)))
            state.set_pressure_oil(Layer.pressure_oil_init * np.ones(Layer.dim, 1))
            state.set_pressure_cap()
            state.set_pressure_water()

    def initialize_cells(self):
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    self.initialize_cell(self.container[k, i, j])

    def equate_cell_states(self, to_previous=False):
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    cell = self.container[k, i, j]
                    if not to_previous:
                        cell.get_cell_state_n().set_equals_to(cell.get_cell_state_n_plus())
                    else:
                        cell.get_cell_state_n_plus().set_equals_to(cell.get_cell_state_n())
