import numpy as np

from .Layer import Layer
from Enums import Components


class ThreeDimOilWaterImpes:
    def __init__(self):
        self.tau_default = 86400
        self.tau = self.tau_default
        self.delta_0 = 1000  # Начальное приближение, на которое отличаются давления
        self.delta_max = 1

    def generate_delta_k(self):
        delta_k = self.delta_0 * np.ones((Layer.N_z, Layer.N_x, Layer.N_z))
        return delta_k

    def count_norm(self, delta_k):
        return np.amax(delta_k)

    def check_norm(self, delta_k):
        isMore = True
        norm = self.check_norm(delta_k)
        if norm <= self.delta_max:
            isMore = False
        return isMore

    @staticmethod
    def recount_properties(cell_container):
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    cell = cell_container.get_cell(k, i, j)
                    cell_state_n_plus = cell.get_cell_state_n_plus()
                    for component in Layer.components:
                        component_index = Layer.components.index(component)
                        ro = component.count_ro(cell_state_n_plus.get_components_pressure()[component_index])
                        k_r = component.count_k_r(cell_state_n_plus.get_components_saturation()[component_index])
                        cell_state_n_plus.set_ro(ro)
                        cell_state_n_plus.set_k_r(k_r)
                    fi = Layer.count_fi(cell_state_n_plus.get_pressure_oil())
                    cell_state_n_plus.set_fi(fi)

    def count_flows(self, flows):
        for k in range(Layer.N_z - 1):
            for i in range(Layer.N_x - 1):
                for j in range(Layer.N_z - 1):
                    flow = flows[k, i, j]
                    flow.count_flows()

    # TODO: реализовать функции подсчета коэффициентов
    def count_a(self):
        pass

    def count_b(self):
        pass

    def count_c(self):
        pass

    def count_d(self):
        pass

    def count_e(self):
        pass

    def count_f(self):
        pass

    def count_g(self):
        pass

    def count_h(self):
        pass



    def generate_matrix(self):
        # TODO: Собственно сгенерировать всю матрицу. Ну которая семидиагональная ага да
        pass

    def solve_slau(self):
        # TODO: Реши систему. Тащемта можешь использовать стандартный решатель
        pass

    def update_pressure(self, cell_container, delta_k):
        # TODO: Обнови давление ага да
        pass

    def update_saturation(self, cell_container, flows):
        # TODO: здесь обнови насыщенность
        pass

    def update_pressure_cap(self, cell_container):
        # TODO: КАПИЛЯРКА АГА (обнови) Перестань общаться с собой через комментарии к коду!!
        pass

    def check_pressure_convergence(self, cell_container, delta_k):
        # TODO: Проверить сходимость по давлению
        pass