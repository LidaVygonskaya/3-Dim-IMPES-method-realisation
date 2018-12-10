import numpy as np

from Layer import Layer
from Cell import Cell
from Enums import Components


class Flow:
    def __init__(self, left_cell, right_cell):
        self.t_oil_water = np.zeros((Layer.dim, Layer.components_count))
        self.left_cell = left_cell
        self.right_cell = right_cell

    @staticmethod
    def initialize_flow_array(cell_container):
        """
        Инициализирует потоки. Проставляет им левые и правые ячейки. В случае граничных ячеек потоки равен нулю.
        В случае граничных ячеек пока что у потока стоит правая клетка равная левой клетке.
        Размер матрицы потоков на один меньше с каждой стороны чем основного массива
        :param cell_container: контрейнер с клетками. Объект класса CellContainer
        :return: массив потоков. Массив объектов типа Flow
        """
        flow_array = np.zeros((Layer.N_z, Layer.N_x, Layer.N_y), dtype=Flow)
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    cell = cell_container.get_cell(k, i, j)
                    left_cell = np.array([cell,
                                          cell,
                                          cell])

                    right_cell_x = None
                    right_cell_y = None
                    right_cell_z = None

                    if i != Layer.N_x - 1:
                        right_cell_x = cell_container.get_cell(k, i + 1, j)

                    if j != Layer.N_y - 1:
                        right_cell_y = cell_container.get_cell(k, i, j + 1)

                    if k != Layer.N_z - 1:
                        right_cell_z = cell_container.get_cell(k + 1, i, j)

                    right_cell = np.array([right_cell_x,
                                          right_cell_y,
                                          right_cell_z])

                    flow_array[k, i, j] = Flow(left_cell, right_cell)

        return flow_array

    def get_max_pressure_cell(self, dimIndex, componentIndex):
        if self.right_cell[dimIndex] is None:
            return self.left_cell[dimIndex]
        left_cell_pressure = self.left_cell[dimIndex].get_cell_state_n_plus().get_components_pressure()
        right_cell_pressure = self.right_cell[dimIndex].get_cell_state_n_plus().get_components_pressure()
        if left_cell_pressure[componentIndex] >= right_cell_pressure[componentIndex]:
            return self.left_cell[dimIndex]
        else:
            return self.right_cell[dimIndex]


    # TODO: зови Яна. Одной лень
    def count_flow(self):
        for componentIndex in range(Layer.components_count):
            t_component_table = np.zeros(Layer.dim)  # Столбец размером Layer.dim. Будем заполнять матрицу по столбцам для каждой компоненты
            for dimIndex in range(Layer.dim):
                cell = self.get_max_pressure_cell(dimIndex, componentIndex)
                if cell is None:
                    t_component_element = 0
                else:
                    cell_state_n_plus = cell.get_cell_state_n_plus()
                    # TODO: По сути k - это вектор, поэтому нужно будет что-то с ним сделать
                    t_component_element = (cell.get_k(dimIndex) * cell_state_n_plus.get_components_k_r()[componentIndex] /
                                       cell.get_mu_oil_water()[componentIndex]) * (1 / Layer.h) ** 2.0 * \
                                       cell_state_n_plus.get_components_ro()[componentIndex]
                t_component_table[dimIndex] = t_component_element
            self.t_oil_water[:, componentIndex] = t_component_table

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