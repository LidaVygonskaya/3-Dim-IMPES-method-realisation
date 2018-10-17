import numpy as np
from .Layer import Layer


class Flow:
    def __init__(self, left_cell, right_cell):
        self.left_cell = left_cell
        self.right_cell = right_cell
        self.t_oil_water = np.zeros((3, 2))

    @staticmethod
    def initialize_flow_array():
        # TODO: создать массив потоков. Он должен быть размером N_x-1 * N_y-1 * N_z-1
        # TODO: Когда будешь пропихивать ему две клетки посмотри индексы. Просто афигеть как важно
        pass

    def get_right_cell(self):
        return self.right_cell

    def get_left_cell(self):
        return self.left_cell

    def get_max_pressure_cell(self, index):
        # TODO: получить клетку вверх по потоку
        pass

    def count_flows(self):
        # TODO: поссчитать все потоки
        for i in range(Layer.components_count):
            cell = self.get_max_pressure_cell(i)
            cell_state_n_plus = cell.get_cell_state_n_plus()
            t_component = (cell.get_k() * cell_state_n_plus.get_components_k_r()[i] / cell.get_mu_oil_water()[i]) \
                  * (1 / cell.layer.h) ** 2.0 * cell_state_n_plus.get_components_ro()[i]
            self.t_oil_water[:, i] = t_component