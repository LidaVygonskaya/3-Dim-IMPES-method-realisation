import numpy as np
from scipy.sparse import bsr_matrix
from Layer import Layer


class SolverSlau:
    def __init__(self):
        self.N_x = Layer.N_x
        self.N_y = Layer.N_y
        self.N_z = Layer.N_z
        self.dim = Layer.N_x * Layer.N_y * Layer.N_z

        #self.coefficient_matrix = np.array((self.dim, self.dim))
        #self.coefficient_matrix = bsr_matrix((self.dim, self.dim))
        # TODO: sparse diags посмотри пример придется все переделывать
        self.nevyaz_vector = np.zeros(self.dim)
        self.result_vector = []

    def add_coefficient(self, i, j, coefficient):
        self.coefficient_matrix[i, j] += coefficient

    def set_coefficient(self, i, j, value):
        self.coefficient_matrix[i, j] = value

    def set_zero(self):
        self.coefficient_matrix = np.zeros((self.dim, self.dim))
        self.nevyaz_vector = np.zeros(self.dim)

    def add_nevyaz(self, i, coefficient):
        self.nevyaz_vector[i] += coefficient

    def get_shift_N_x(self):
        return self.N_x

    def get_shift_N_xy(self):
        return self.N_x * self.N_y

    def set_matrix_coefficients(self, i, j, coefficient):
        if self.dim > i >= 0:
            self.add_coefficient(i, i, coefficient)  # Добавляем этот же коэффициент на место a
            # TODO: организовать добавление в невязку этого коэффициента
            # TODO: смотри аналогичное добавление этого в предыдущем коде
            if self.dim > j >= 0:
                self.add_coefficient(i, j, coefficient)  # Добавляем коэффициент в матрицу


