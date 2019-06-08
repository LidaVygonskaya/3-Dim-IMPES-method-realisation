import numpy as np

from Layer import Layer
from Enums import Components, FlowComponents


class Flow:
    def __init__(self, left_cell=None, right_cell=None):
        self.t_oil_water = np.zeros(Layer.components_count)
        self.left_cell = left_cell
        self.right_cell = right_cell

    @staticmethod
    def initialize_flow(cell_container):
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    cell = cell_container.get_cell(k, i, j)

                    flow_array_x = cell.get_flow('x')
                    flow_array_y = cell.get_flow('y')
                    flow_array_z = cell.get_flow('z')

                    for component_index in range(Layer.components_count):
                        #========================= Для x ===============================================================
                        if i == Layer.N_x - 1:
                            flow_array_x[FlowComponents.minus.value].left_cell = cell_container.get_cell(k, i - 1, j)
                            flow_array_x[FlowComponents.minus.value].right_cell = cell

                            flow_array_x[FlowComponents.plus.value].left_cell = cell
                            flow_array_x[FlowComponents.plus.value].right_cell = None

                        elif i == 0:
                            flow_array_x[FlowComponents.minus.value].left_cell = None
                            flow_array_x[FlowComponents.minus.value].right_cell = cell

                            flow_array_x[FlowComponents.plus.value].left_cell = cell
                            flow_array_x[FlowComponents.plus.value].right_cell = cell_container.get_cell(k, i + 1, j)

                        else:
                            flow_array_x[FlowComponents.minus.value].left_cell = cell_container.get_cell(k, i - 1, j)
                            flow_array_x[FlowComponents.minus.value].right_cell = cell

                            flow_array_x[FlowComponents.plus.value].left_cell = cell
                            flow_array_x[FlowComponents.plus.value].right_cell = cell_container.get_cell(k, i + 1, j)

                        #============================= Для y ===========================================================
                        if j == Layer.N_y - 1:
                            flow_array_y[FlowComponents.minus.value].left_cell = cell_container.get_cell(k, i, j - 1)
                            flow_array_y[FlowComponents.minus.value].right_cell = cell

                            flow_array_y[FlowComponents.plus.value].left_cell = cell
                            flow_array_y[FlowComponents.plus.value].right_cell = None

                        elif j == 0:
                            flow_array_y[FlowComponents.minus.value].left_cell = None
                            flow_array_y[FlowComponents.minus.value].right_cell = cell

                            flow_array_y[FlowComponents.plus.value].left_cell = cell
                            flow_array_y[FlowComponents.plus.value].right_cell = cell_container.get_cell(k, i, j + 1)
                        else:
                            flow_array_y[FlowComponents.minus.value].left_cell = cell_container.get_cell(k, i, j - 1)
                            flow_array_y[FlowComponents.minus.value].right_cell = cell

                            flow_array_y[FlowComponents.plus.value].left_cell = cell
                            flow_array_y[FlowComponents.plus.value].right_cell = cell_container.get_cell(k, i, j + 1)

                        #============================= Для z ===========================================================
                        if k == Layer.N_z - 1:
                            flow_array_z[FlowComponents.minus.value].left_cell = cell_container.get_cell(k - 1, i, j)
                            flow_array_z[FlowComponents.minus.value].right_cell = cell

                            flow_array_z[FlowComponents.plus.value].left_cell = cell
                            flow_array_z[FlowComponents.plus.value].right_cell = None
                        elif k == 0:
                            flow_array_z[FlowComponents.minus.value].left_cell = None
                            flow_array_z[FlowComponents.minus.value].right_cell = cell

                            flow_array_z[FlowComponents.plus.value].left_cell = cell
                            flow_array_z[FlowComponents.plus.value].right_cell = cell_container.get_cell(k + 1, i, j)
                        else:
                            flow_array_z[FlowComponents.minus.value].left_cell = cell_container.get_cell(k - 1, i, j)
                            flow_array_z[FlowComponents.minus.value].right_cell = cell

                            flow_array_z[FlowComponents.plus.value].left_cell = cell
                            flow_array_z[FlowComponents.plus.value].right_cell = cell_container.get_cell(k + 1, i, j)
                        #===============================================================================================

    def get_max_pressure_cell(self, component_index):
        left_cell = self.left_cell
        right_cell = self.right_cell

        if left_cell is None:
            return right_cell
        elif right_cell is None:
            return left_cell
        else:
            left_cell_pressure = self.left_cell.get_cell_state_n_plus().get_components_pressure()
            right_cell_pressure = self.right_cell.get_cell_state_n_plus().get_components_pressure()
            if left_cell_pressure[component_index] >= right_cell_pressure[component_index]:
                return left_cell
            else:
                return right_cell

    def count_flow(self, direction):
        for component_index in range(Layer.components_count):
            t_component = 0
            cell = self.get_max_pressure_cell(component_index)  # Max pressure cell
            cell_state_n_plus = cell.get_cell_state_n_plus()
            if self.left_cell is not None or self.right_cell is not None:
                t_component = (cell.get_k(direction) * cell_state_n_plus.get_components_k_r()[component_index]
                               / cell.get_mu_oil_water()[component_index]) * Layer.V_ijk\
                              * (1 / Layer.get_h(direction)) ** 2.0 * cell_state_n_plus.get_components_ro()[component_index]
            self.t_oil_water[component_index] = t_component

    def get_left_cell(self):
        return self.left_cell

    def get_right_cell(self):
        return self.right_cell

    def get_water_flow_vector(self):
        return self.t_oil_water[:, Components.WATER.value]

    def get_oil_flow_vector(self):
        return self.t_oil_water[:, Components.OIL.value]

    def get_water_flow(self):
        return self.t_oil_water[Components.WATER.value]

    def get_oil_flow(self):
        return self.t_oil_water[Components.OIL.value]
