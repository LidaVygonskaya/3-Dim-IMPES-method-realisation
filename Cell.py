from CellState import CellState
from Layer import Layer


class Cell:
    def __init__(self, eq_index, is_boundary_x=False, is_boundary_y=False):
        self.is_boundary_y = is_boundary_y
        self.is_boundary_x = is_boundary_x
        self.cell_states = [CellState(), CellState()]#n, n + 1 layers
        self.eq_index = eq_index

    def is_n_plus_state(self, state):
        if self.cell_states.index(state) == 1:
            return True
        else:
            return False

    def is_boundary_cell_y(self):
        """
        Проверяет является ли клетка граничной
        :return: True - если клетка граничная. False - если нет
        """
        return self.is_boundary_y

    def is_boundary_cell_x(self):
        """
        Проверяет является ли клетка граничной
        :return: True - если клетка граничная. False - если нет
        """
        return self.is_boundary_x

    # TODO: rewrite all functions
    #Cell constant parameters
    def get_eq_index(self):
        return self.eq_index

    def get_k(self, dimIndex):
        """
        Возвращает проницаемость определенной компоненты x, y или z
        :param dimIndex:
        :return: k
        """
        return Layer.k_array[dimIndex]

    def get_c_f_water(self):
        pass

    def get_c_f_oil(self):
        pass

    def get_mu_oil(self):
        pass

    def get_mu_water(self):
        pass

    def get_mu_oil_water(self):
        return Layer.mu_oil_water

    def get_cell_state_n(self):
        return self.cell_states[0]

    def get_cell_state_n_plus(self):
        return self.cell_states[1]

    def get_cell_states(self):
        return self.cell_states