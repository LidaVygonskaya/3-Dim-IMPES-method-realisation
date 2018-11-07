from CellState import CellState


class Cell:
    def __init__(self, eq_index, is_boundary=False):
        self.is_boundary = is_boundary
        self.cell_states = [CellState(), CellState()]#n, n + 1 layers
        self.eq_index = eq_index

    def is_boundary_cell(self):
        """
        Проверяет является ли клетка граничной
        :return: True - если клетка граничная. False - если нет
        """
        return self.is_boundary

    # TODO: rewrite all functions
    #Cell constant parameters
    def get_eq_index(self):
        return self.eq_index

    def get_k(self):
        return self.layer.k

    def get_c_f_water(self):
        return self.layer.c_f_water

    def get_c_f_oil(self):
        return self.layer.c_f_oil

    def get_mu_oil(self):
        return self.layer.mu_oil

    def get_mu_water(self):
        return self.layer.mu_water

    def get_mu_oil_water(self):
        return self.layer.mu_oil_water

    def get_cell_state_n(self):
        return self.cell_states[0]

    def get_cell_state_n_plus(self):
        return self.cell_states[1]

    def get_cell_states(self):
        return self.cell_states