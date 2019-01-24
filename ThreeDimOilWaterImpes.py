import numpy as np
from scipy.sparse import diags

from Layer import Layer
from Enums import FlowComponents


class ThreeDimOilWaterImpes:
    def __init__(self, solver_slau):
        self.tau_default = 86400
        self.tau = self.tau_default
        self.time_max = self.tau_default * 365
        self.delta_0 = 1000  # Начальное приближение, на которое отличаются давления
        self.delta_max = 10 ** (-3)
        self.solver_slau = solver_slau

    def generate_delta_k(self):
        delta_k = self.delta_0 * np.ones(Layer.N_z * Layer.N_x * Layer.N_z)
        return delta_k

    def count_norm(self, delta_k):
        return np.amax(np.absolute(delta_k))

    def check_norm(self, delta_k):
        isMore = True
        norm = self.count_norm(delta_k)
        print(norm)
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
    @staticmethod
    def count_cells_flows(cell_container):
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    cell = cell_container.get_cell(k, i, j)
                    flow_array_x = cell.get_flow('x')
                    flow_array_y = cell.get_flow('y')
                    flow_array_z = cell.get_flow('z')

                    for flow in flow_array_x:
                        flow.count_flow('x')

                    for flow in flow_array_y:
                        flow.count_flow('y')

                    for flow in flow_array_z:
                        flow.count_flow('z')

    def count_c1_p(self, cell):
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        return ((state_n.get_fi() * state_n.get_s_water() * Layer.water.ro_water_0 * Layer.water.c_f_water) + (
            state_n.get_s_water() * state_n_plus.get_ro_water() * Layer.c_r * Layer.fi_0)) / self.tau

    def count_c2_p(self, cell):
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        return ((1.0 - state_n.get_s_water()) * (
            state_n.get_fi() * Layer.oil.ro_oil_0 * Layer.oil.c_f_oil + Layer.c_r * Layer.fi_0 * state_n_plus.get_ro_oil())) / self.tau

    # TODO: реализовать функции подсчета коэффициентов
    def count_a(self, cell):
        # Коэффициент в матрицу потоков берется просто как разность всех коэффициентов
        # Здесь считается только то, что пойдет в матрицу потоков не от други кожффициентов
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        #p_cap_graph = cell.layer.count_pcap_graph()
        #p_cap_der = cell.layer.count_p_cap_graph_der(p_cap_graph, state_n.get_s_water())
        #ro_der = cell.layer.ro_water_0 * cell.layer.c_f_water
        A = 1 # TODO: свериться. На самом деле помоему коэффициент не совесем такой
        # Вот такой коэффициент
        coeff = A * self.count_c1_p(cell) + self.count_c2_p(cell)
        # TODO: set coefficients a_d similar as in filtration и пререписать c1_p и c2_p
        return coeff

    def count_coefficient(self, flow, dim_string):
        # TODO: реализовать одну функцию для подсчета всех коэффициентов
        pass

    def count_b(self, flow):
        # Плюсовый коэффициент
        # Стоит в матрице на (x, y + Nx * Ny). P(i+1, j, k)
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()
        if right_cell is None:
            return 0
        pressure_left = left_cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_right = right_cell.get_cell_state_n_plus().get_pressure_oil()
        A = 1
        b = A * flow.get_water_flow() + flow.get_oil_flow()
        r_ost = -A * flow.get_water_flow() * (
            right_cell.get_cell_state_n_plus().get_pressure_cap() - left_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(right_cell.get_eq_index(), -r_ost - b * (pressure_left - pressure_right))
        return b

    def count_c(self, flow):
        # Стоит в матрице на (x + Nx * Ny, y). P(i-1, j, k)
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()
        if left_cell is None:
            return 0
        pressure_left = left_cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_right = right_cell.get_cell_state_n_plus().get_pressure_oil()
        A = 1
        c = A * flow.get_water_flow() + flow.get_oil_flow()
        r_ost = -A * flow.get_water_flow() * (
            left_cell.get_cell_state_n_plus().get_pressure_cap() - right_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(left_cell.get_eq_index(), -r_ost - c * (pressure_right - pressure_left))
        return c

    def count_d(self, flow):
        # Стоит в матрице на (x + Nx * Ny, y). P(i, j + 1, k)
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()
        if right_cell is None:
            return 0
        pressure_left = left_cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_right = right_cell.get_cell_state_n_plus().get_pressure_oil()
        A = 1
        d = A * flow.get_water_flow() + flow.get_oil_flow()
        r_ost = -A * flow.get_water_flow() * (
            right_cell.get_cell_state_n_plus().get_pressure_cap() - left_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(right_cell.get_eq_index(), -r_ost - d * (pressure_left - pressure_right))
        return d

    def count_e(self, flow):
        # Стоит в матрице на (x + Nx * Ny, y). P(i, j - 1, k)
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()
        A = 1
        e = A * flow.get_water_flow() + flow.get_oil_flow()
        if left_cell is None:
            return 0
        pressure_left = left_cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_right = right_cell.get_cell_state_n_plus().get_pressure_oil()
        r_ost = -A * flow.get_water_flow() * (
            left_cell.get_cell_state_n_plus().get_pressure_cap() - right_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(left_cell.get_eq_index(), -r_ost - e * (pressure_right - pressure_left))
        return e

    def count_g(self, flow):
        # Стоит в матрице на (x + Nx * Ny, y). P(i, j, k + 1)
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()
        if right_cell is None:
            return 0
        pressure_left = left_cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_right = right_cell.get_cell_state_n_plus().get_pressure_oil()
        A = 1
        g = A * flow.get_water_flow() + flow.get_oil_flow()
        r_ost = -A * flow.get_water_flow() * (
            right_cell.get_cell_state_n_plus().get_pressure_cap() - left_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(right_cell.get_eq_index(), -r_ost - g * (pressure_left - pressure_right))
        return g

    def count_f(self, flow):
        # Стоит в матрице на (x + Nx * Ny, y). P(i, j, k - 1)
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()
        if left_cell is None:
            return 0
        A = 1
        f = A * flow.get_water_flow() + flow.get_oil_flow()
        pressure_left = left_cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_right = right_cell.get_cell_state_n_plus().get_pressure_oil()
        r_ost = -A * flow.get_water_flow() * (
            left_cell.get_cell_state_n_plus().get_pressure_cap() - right_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(left_cell.get_eq_index(), -r_ost - f * (pressure_right - pressure_left))
        return f

    def generate_matrix(self, cell_container):
        # TODO: Собственно сгенерировать всю матрицу. Ну которая семидиагональная ага да
        # TODO: Посчитать коэффициент a
        matrix_size = Layer.N_x * Layer.N_y * Layer.N_z
        b = np.empty((0, matrix_size))  # Плюсовый
        c = np.empty((0, matrix_size))  # Минусовый
        d = np.empty((0, matrix_size))  # Плюсовый
        e = np.empty((0, matrix_size))  # Минусовый
        g = np.empty((0, matrix_size))  # Плюсовый
        f = np.empty((0, matrix_size))  # Минусовый
        coeff_a_d = np.empty((0, matrix_size))
        coeff_a_d_nevyaz = np.empty((0, matrix_size))

        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    #if i == 3 and k == 3 and j == 3:
                    #    print("be")
                    cell = cell_container.get_cell(k, i, j)
                    pressure_n = cell.get_cell_state_n().get_pressure_oil()
                    pressure_n_plus = cell.get_cell_state_n_plus().get_pressure_oil()

                    b = np.append(b, self.count_b(cell.get_flow_coefficient('x', FlowComponents.plus.value)))  # Откусываем большой сдвиг с конца
                    c = np.append(c, self.count_c(cell.get_flow_coefficient('x', FlowComponents.minus.value)))  # Откусываем большой сдвиг с начала
                    d = np.append(d, self.count_d(cell.get_flow_coefficient('y', FlowComponents.plus.value)))  # Откусываем малый сдвиг с конца
                    e = np.append(e, self.count_e(cell.get_flow_coefficient('y', FlowComponents.minus.value)))  # Откусываем малый сдвиг с начала
                    g = np.append(g, self.count_g(cell.get_flow_coefficient('z', FlowComponents.plus.value)))  # Остается неизменным сдвиг + 1
                    f = np.append(f, self.count_f(cell.get_flow_coefficient('z', FlowComponents.minus.value)))  # Остается неизменным сдвиг -1

                    # Добавляем в середину
                    a_d = self.count_a(cell)
                    coeff_a_d = np.append(coeff_a_d, a_d)
                    coeff_a_d_nevyaz = np.append(coeff_a_d_nevyaz, a_d * (pressure_n_plus - pressure_n))

        #TODO:добавляем в невязку
        # Невязка a_d
        self.solver_slau.add_vector_to_nevyaz(np.transpose([coeff_a_d_nevyaz]))
        # Невязка от ф
        a = - b - c - d - e - f - g - coeff_a_d
        # TODO: откусить
        shift_Nx = self.solver_slau.get_shift_N_x()
        shift_Nxy = self.solver_slau.get_shift_N_xy()

        f = f[(shift_Nxy):]
        g = g[:-(shift_Nxy)]

        e = e[1:]
        d = d[:-1]

        c = c[(shift_Nx):]
        b = b[:(-shift_Nx)]

        diagonals = [a, e, d, c, b, f, g]
        shifts = [0, -1, 1, -shift_Nx, shift_Nx, -shift_Nxy, shift_Nxy]
        self.solver_slau.coefficient_matrix = diags(diagonals, shifts)
        # Для смотрения
        smotrenie = self.solver_slau.coefficient_matrix.toarray()
        print("be")

    def solve_slau(self):
       self.solver_slau.solve_slau()

    def update_pressure(self, cell_container, delta_k):
        eq_index = 0
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    cell = cell_container.get_cell(k, i, j)
                    state_n_plus = cell.get_cell_state_n_plus()
                    state_n_plus.set_pressure_oil(state_n_plus.get_pressure_oil() + delta_k[eq_index])
                    state_n_plus.set_pressure_water(state_n_plus.get_pressure_water() - state_n_plus.get_pressure_cap())
                    eq_index += 1

    def update_saturation(self, cell_container):
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    # TODO: получать разные клетки по всем осям
                    cell = cell_container.get_cell(k, i, j)
                    flow_x_plus = cell.get_flow_coefficient('x', FlowComponents.plus.value)
                    flow_x_minus = cell.get_flow_coefficient('x', FlowComponents.minus.value)

                    flow_y_plus = cell.get_flow_coefficient('y', FlowComponents.plus.value)
                    flow_y_minus = cell.get_flow_coefficient('y', FlowComponents.minus.value)

                    flow_z_plus = cell.get_flow_coefficient('z', FlowComponents.plus.value)
                    flow_z_minus = cell.get_flow_coefficient('z', FlowComponents.minus.value)

                    cell_x_plus = cell_container.get_cell(k, i + 1, j)
                    cell_x_minus = cell_container.get_cell(k, i - 1, j)

                    cell_y_plus = cell_container.get_cell(k, i, j + 1)
                    cell_y_minus = cell_container.get_cell(k, i, j - 1)

                    cell_z_plus = cell_container.get_cell(k + 1, i, j)
                    cell_z_minus = cell_container.get_cell(k - 1, i, j)

                    state_n = cell.get_cell_state_n()
                    state_n_plus = cell.get_cell_state_n_plus()
                    #================ water ================
                    # TODO: Запихать эти штуки класс для воды или нефти
                    ro_der = Layer.water.ro_water_0 * Layer.water.c_f_water
                    fi_der = Layer.fi_0 * Layer.c_r

                    p_cap_der = Layer.count_p_cap_graph_der(state_n.get_s_water())

                    d_11 = (1.0 / self.tau) * state_n.get_s_water() * (
                    state_n.get_fi() * ro_der + fi_der * state_n_plus.get_ro_water())
                    d_12 = (1.0 / self.tau) * (
                    state_n_plus.get_fi() * state_n_plus.get_ro_water() - ro_der * state_n.get_s_water() * state_n.get_ro_water() * p_cap_der)

                    coeff_2 = 0

                    coeff_3 = -(d_11 / d_12) * (state_n_plus.get_pressure_oil() - state_n.get_pressure_oil())

                    s_water_new = state_n.get_s_water() + (1.0 / d_12) * coeff_2 + coeff_3

                    state_n_plus.set_s_water(s_water_new)
                    state_n_plus.set_s_oil(1.0 - s_water_new)

    # TODO: ты чето тут сломала посмотри
    def update_pressure_cap(self, cell_container):
        for k in range(1, Layer.N_z - 1):
            for i in range(1, Layer.N_x - 1):
                for j in range(1, Layer.N_y - 1):
                    directions = ['x', 'y', 'z']
                    cell = cell_container.get_cell(k, i, j)

                    s_water = cell.get_cell_state_n_plus().get_s_water()
                    p_cap = Layer.count_pressure_cap(s_water)
                    cell.get_cell_state_n_plus().set_pressure_cap(p_cap)

    def check_pressure_convergence(self, cell_container, delta_k):
        # TODO: Проверить сходимость по давлению не факт что оно нужно
        pass