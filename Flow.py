import numpy as np

from Layer import Layer
from Cell import Cell
from Enums import Components


class Flow:
    def __init__(self, left_cell, right_cell):
        self.t_oil_water = np.zeros((Layer.dim, Layer.components_count))
        self.left_cell = left_cell#np.zeros((Layer.dim, Layer.components_count), dtype=Cell)
        self.right_cell = right_cell#np.zeros((Layer.dim, Layer.components_count), dtype=Cell)

    @staticmethod
    def initialize_flow_array(cell_container):
        """
        Инициализирует потоки. Проставляет им левые и правые ячейки. В случае граничных ячеек потоки равен нулю.
        В случае граничных ячеек пока что у потока стоит правая клетка равная левой клетке
        :param cell_container: контрейнер с клетками. Объект класса CellContainer
        :return: массив потоков. Массив объектов типа Flow
        """
        flow_array = np.zeros((Layer.N_z, Layer.N_x, Layer.N_y), dtype=Flow)
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    cell = cell_container.get_cell(k, i, j)
                    if cell.is_boundary_cell():
                        flow_array[k, i, j] = Flow(cell, cell)
                    else:
                        left_cell = np.array([cell_container.get_cell(k, i, j),
                                              cell_container.get_cell(k, i, j),
                                              cell_container.get_cell(k, i, j)])

                        right_cell = np.array([cell_container.get_cell(k, i + 1, j),
                                              cell_container.get_cell(k, i, j + 1),
                                              cell_container.get_cell(k + 1, i, j)])

                        flow_array[k, i, j] = Flow(left_cell, right_cell)

        return flow_array

    def get_max_pressure_cell(self, index):
        # TODO: получить клетку вверх по потоку
        pass

    def count_flow(self):
        # TODO: поссчитать один поток
        for i in range(Layer.components_count):
            cell = self.get_max_pressure_cell(i)
            cell_state_n_plus = cell.get_cell_state_n_plus()
            t_component = (cell.get_k() * cell_state_n_plus.get_components_k_r()[i] / cell.get_mu_oil_water()[i]) \
                  * (1 / cell.layer.h) ** 2.0 * cell_state_n_plus.get_components_ro()[i]
            self.t_oil_water[:, i] = t_component

    def get_water_flow_vector(self):
        return self.t_oil_water[:, Components.WATER.value]

    def get_oil_flow_vector(self):
        return self.t_oil_water[:, Components.OIL.value]

    def get_right_cell(self, dim):
        if dim == 'x':
            return self.right_cell[0]
        elif dim == 'y':
            return self.right_cell[1]
        elif dim == 'z':
            return self.right_cell[2]

    def get_left_cell(self, dim):
        if dim == 'x':
            return self.left_cell[0]
        elif dim == 'y':
            return self.left_cell[1]
        elif dim == 'z':
            return self.left_cell[2]

    def get_water_flow(self, dim):
        if dim == 'x':
            return self.t_oil_water[0, Components.WATER.value]

        elif dim == 'y':
            return self.t_oil_water[1, Components.WATER.value]

        elif dim == 'z':
            return self.t_oil_water[2, Components.WATER.value]

    def get_oil_flow(self, dim):
        if dim == 'x':
            return self.t_oil_water[0, Components.OIL.value]

        elif dim == 'y':
            return self.t_oil_water[1, Components.OIL.value]

        elif dim == 'z':
            return self.t_oil_water[2, Components.OIL.value]