import numpy as np
from scipy.sparse import diags

from Layer import Layer


class ThreeDimOilWaterImpes:
    def __init__(self, solver_slau):
        self.tau_default = 86400
        self.tau = self.tau_default
        self.delta_0 = 1000  # Начальное приближение, на которое отличаются давления
        self.delta_max = 1
        self.solver_slau = solver_slau

    def generate_delta_k(self):
        delta_k = self.delta_0 * np.ones((Layer.N_z, Layer.N_x, Layer.N_z))
        return delta_k

    def count_norm(self, delta_k):
        return np.amax(delta_k)

    def check_norm(self, delta_k):
        isMore = True
        norm = self.count_norm(delta_k)
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
                        cell_state_n_plus.set_ro(ro, component_index)
                        cell_state_n_plus.set_k_r(k_r, component_index)
                    fi = Layer.count_fi(cell_state_n_plus.get_pressure_oil())
                    cell_state_n_plus.set_fi(fi)

    @staticmethod
    def count_flows(flows):
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_z):
                    flow = flows[k, i, j]
                    flow.count_flow()

    # TODO: реализовать функции подсчета коэффициентов
    def count_a(self, cell):
        # Коэффициент в матрицу потоков берется просто как разность всех коэффициентов
        # Здесь считается только то, что пойдет в матрицу потоков не от други кожффициентов
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        p_cap_graph = cell.layer.count_pcap_graph()
        p_cap_der = cell.layer.count_p_cap_graph_der(p_cap_graph, state_n.get_s_water())
        ro_der = cell.layer.ro_water_0 * cell.layer.c_f_water
        A = 1 # TODO: свериться. На самом деле помоему коэффициент не совесем такой
        # Вот такой коэффициент
        coeff = state_n_plus.get_pressure_oil(), state_n.get_pressure_oil(), A * state_n_plus.get_c1_p() + state_n_plus.get_c2_p()
        # TODO: set coefficients a_d similar as in filtration
        return coeff

    def count_b(self, flow):
        # TODO: добавлять r_ost в невязку, и коэффициент домноженный на давление
        # Стоит в матрице на (x, y + Nx * Ny). P(i+1, j, k)
        left_cell = flow.get_left_cell('x')
        right_cell = flow.get_right_cell('x')
        A = 1
        b = A * flow.get_water_flow('x') + flow.get_oil_flow('x')
        r_ost = -A * flow.get_water_flow('x') * (
            right_cell.get_cell_state_n_plus().get_pressure_cap() - left_cell.get_cell_state_n_plus().get_pressure_cap())
        return b

    def count_c(self, flow):
        # Стоит в матрице на (x + Nx * Ny, y). P(i-1, j, k)
        left_cell = flow.get_left_cell('x')
        right_cell = flow.get_right_cell('x')
        A = 1
        c = A * flow.get_water_flow('x') + flow.get_oil_flow('x')
        r_ost = -A * flow.get_water_flow('x') * (
            left_cell.get_cell_state_n_plus().get_pressure_cap() - right_cell.get_cell_state_n_plus().get_pressure_cap())
        return c

    def count_d(self, flow):
        # Стоит в матрице на (x + Nx * Ny, y). P(i, j + 1, k)
        left_cell = flow.get_left_cell('y')
        right_cell = flow.get_right_cell('y')
        A = 1
        d = A * flow.get_water_flow('y') + flow.get_oil_flow('y')
        r_ost = -A * flow.get_water_flow('y') * (
            right_cell.get_cell_state_n_plus().get_pressure_cap() - left_cell.get_cell_state_n_plus().get_pressure_cap())
        return d

    def count_e(self, flow):
        # Стоит в матрице на (x + Nx * Ny, y). P(i, j - 1, k)
        left_cell = flow.get_left_cell('y')
        right_cell = flow.get_right_cell('y')
        A = 1
        e = A * flow.get_water_flow('y') + flow.get_oil_flow('y')
        r_ost = -A * flow.get_water_flow('y') * (
            left_cell.get_cell_state_n_plus().get_pressure_cap() - right_cell.get_cell_state_n_plus().get_pressure_cap())
        return e

    def count_f(self, flow):
        # Стоит в матрице на (x + Nx * Ny, y). P(i, j, k - 1)
        left_cell = flow.get_left_cell('z')
        right_cell = flow.get_right_cell('z')
        A = 1
        f = A * flow.get_water_flow('z') + flow.get_oil_flow('z')
        r_ost = -A * flow.get_water_flow('y') * (
            left_cell.get_cell_state_n_plus().get_pressure_cap() - right_cell.get_cell_state_n_plus().get_pressure_cap())
        return f

    def count_g(self, flow):
        # Стоит в матрице на (x + Nx * Ny, y). P(i, j, k + 1)
        left_cell = flow.get_left_cell('z')
        right_cell = flow.get_right_cell('z')
        A = 1
        g = A * flow.get_water_flow('z') + flow.get_oil_flow('z')
        r_ost = -A * flow.get_water_flow('y') * (
            right_cell.get_cell_state_n_plus().get_pressure_cap() - left_cell.get_cell_state_n_plus().get_pressure_cap())
        return g

    def generate_matrix(self, flow_array, cell_container, solver_slau):
        # TODO: Собственно сгенерировать всю матрицу. Ну которая семидиагональная ага да
        # TODO: Посчитать коэффициент a
        matrix_size = Layer.N_x * Layer.N_y * Layer.N_z
        b = np.empty((0, matrix_size))
        c = np.empty((0, matrix_size))
        d = np.empty((0, matrix_size))
        e = np.empty((0, matrix_size))
        f = np.empty((0, matrix_size))
        g = np.empty((0, matrix_size))
        coeff_a_d = np.empty((0, matrix_size))
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    flow = flow_array[k, i, j]
                    cell = cell_container.get_cell(k, i, j)
                    b = np.append(b, self.count_b(flow))  # Откусываем большой сдвиг с конца
                    c = np.append(c, self.count_c(flow))  # Откусываем большой сдвиг с начала
                    d = np.append(d, self.count_d(flow))  # Откусываем малый сдвиг с конца
                    e = np.append(e, self.count_e(flow))  # Откусываем малый сдвиг с начала
                    f = np.append(f, self.count_f(flow))  # Остается неизменным сдвиг -1
                    g = np.append(g, self.count_g(flow))  # Остается неизменным сдвиг + 1
                    coeff_a_d = np.append(coeff_a_d, self.count_a(cell))

        a = - b - c - d - e - f - g - coeff_a_d

        # TODO: откусить
        shift_Nx = solver_slau.get_shift_N_x()
        shift_Nxy = solver_slau.get_shift_N_xy()
        b = b[:-(shift_Nxy)]  # Откусываем большой сдвиг с конца
        c = c[(shift_Nxy):]  # Откусываем большой сдвиг с начала
        d = d[:-(shift_Nx)]  # Откусываем малый сдвиг с конца
        e = e[(shift_Nx):]  # Откусываем малый сдвиг с начал
        f = f[1:]
        g = g[:-1]
        diagonals = [a, f, g, e, d, c, b]
        shifts = [0, -1, 1, -shift_Nx, shift_Nx, -shift_Nxy, shift_Nxy]
        solver_slau.coefficient_matrix = diags(diagonals, shifts).toarray()

        print('bebe')




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