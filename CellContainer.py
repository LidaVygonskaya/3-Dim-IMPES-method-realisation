import numpy as np

from Layer import Layer
from Cell import Cell


class CellContainer:
    def __init__(self):
        self.container = np.zeros((Layer.N_z, Layer.N_x, Layer.N_z), dtype=Cell)  # (z, x, y) x = rows, y - columns
        eq_index = 0
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    if self.is_well_cell(k, i, j):
                        self.container[k, i, j] = Cell(eq_index, has_well=True)
                    else:
                        self.container[k, i, j] = Cell(eq_index)
                    eq_index += 1

    def get_cells(self):
        return self.container

    def is_well_cell(self, k, i, j):
        if (k, i, j) in Layer.wells:
            return True
        else:
            return False

    def get_cell(self, k, i, j):
        return self.container[k, i, j]

    @staticmethod
    def initialize_cell(cell):
        for state in cell.get_cell_states():
            state.set_s_water(Layer.s_water_init)
            state.set_s_oil(1.0 - Layer.s_water_init)
            if cell.is_n_plus_state(state):
                state.set_pressure_oil(Layer.pressure_oil_init + Layer.delta_0)
            else:
                state.set_pressure_oil(Layer.pressure_oil_init)
            state.set_pressure_cap(Layer.count_pressure_cap(state.get_s_water()))
            state.set_pressure_water(state.get_pressure_oil() - state.get_pressure_cap())  # Отнимается от каждого элемента матрицы

            for i in range(Layer.components_count):
                component = Layer.components[i]
                component_saturation = state.get_components_saturation()[i]
                component_pressure = state.get_components_pressure()[i]
                state.set_k_r(component.count_k_r(component_saturation), i)
                state.set_ro(component.count_ro(component_pressure), i)

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
