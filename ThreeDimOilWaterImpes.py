import numpy as np
from scipy.sparse import diags

from Layer import Layer
from Enums import FlowComponents, Components


class ThreeDimOilWaterImpes:
    def __init__(self, solver_slau):
        self.tau_default = 86400
        self.tau = self.tau_default
        self.time_max = self.tau_default * 365
        self.delta_0 = 1000  # Начальное приближение, на которое отличаются давления
        self.delta_max = 10 ** (-3)
        self.solver_slau = solver_slau
        matrix_size = Layer.N_x * Layer.N_y * Layer.N_z
        self.b = np.empty((0, matrix_size))  # Плюсовый
        self.c = np.empty((0, matrix_size))  # Минусовый
        self.d = np.empty((0, matrix_size))  # Плюсовый
        self.e = np.empty((0, matrix_size))  # Минусовый
        self.g = np.empty((0, matrix_size))  # Плюсовый
        self.f = np.empty((0, matrix_size))  # Минусовый
        self.well_addition_oil = np.empty((0, matrix_size))
        self.well_addition_water = np.empty((0, matrix_size))
        self.coeff_a_d = np.empty((0, matrix_size))
        self.coeff_a_d_nevyaz = np.empty((0, matrix_size))
        self.coeff_well_nevyaz = np.empty((0, matrix_size))

    def set_zero_coeffs(self):
        matrix_size = Layer.N_x * Layer.N_y * Layer.N_z
        self.b = np.empty((0, matrix_size))  # Плюсовый
        self.c = np.empty((0, matrix_size))  # Минусовый
        self.d = np.empty((0, matrix_size))  # Плюсовый
        self.e = np.empty((0, matrix_size))  # Минусовый
        self.g = np.empty((0, matrix_size))  # Плюсовый
        self.f = np.empty((0, matrix_size))  # Минусовый
        self.well_addition_oil = np.empty((0, matrix_size))
        self.well_addition_water = np.empty((0, matrix_size))
        self.coeff_a_d = np.empty((0, matrix_size))
        self.coeff_a_d_nevyaz = np.empty((0, matrix_size))
        self.coeff_well_nevyaz = np.empty((0, matrix_size))

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
    def recount_property(cell_container, k, i, j):
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
        if cell.has_well:
            cell.well.count_well_index(cell, (k, i, j))

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
                    if cell.has_well:
                        cell.well.count_well_index(cell, (k, i, j))

    @staticmethod
    def count_flows(flows):
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_z):
                    flow = flows[k, i, j]
                    flow.count_flow()

    @staticmethod
    def count_cell_flow(cell_container, k, i, j):
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
            state_n.get_s_water() * state_n_plus.get_ro_water() * Layer.c_r * Layer.fi_0)) * Layer.V_ijk / self.tau

    def count_c2_p(self, cell):
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        return ((1.0 - state_n.get_s_water()) * (
            state_n.get_fi() * Layer.oil.ro_oil_0 * Layer.oil.c_f_oil + Layer.c_r * Layer.fi_0 * state_n_plus.get_ro_oil())) * Layer.V_ijk / self.tau

    def count_A(self, state_n, state_n_plus, ):
        p_cap_der = Layer.count_p_cap_graph_der(state_n.get_s_water())
        ro_der = Layer.water.ro_water_0 * Layer.water.c_f_water
        A = (state_n_plus.get_fi() * state_n_plus.get_ro_oil()) / (state_n_plus.get_fi() * state_n_plus.get_ro_water() - state_n.get_s_water() * state_n.get_fi() * ro_der * p_cap_der)
        return A

    # TODO: реализовать функции подсчета коэффициентов
    def count_a(self, cell):
        # TODO: Поменять коэффициент А везде
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        A = self.count_A(state_n, state_n_plus)
        # Вот такой коэффициент
        coeff = A * self.count_c1_p(cell) + self.count_c2_p(cell)
        # TODO: set coefficients a_d similar as in filtration и пререписать c1_p и c2_p
        return coeff

    def count_well_addition(self, cell, water=False):
        # TODO: посмотреть для нагнетательной скавжины насыщенность и откуда она берется
        # TODO: оно должно браться из класса скважины
        """
        Добавка к коэффициентам от скважинки
        :param cell:
        :return:
        """
        if cell.has_well:
            well_index = cell.well.well_index_oil_water
            if water:
                state_n = cell.get_cell_state_n()
                state_n_plus = cell.get_cell_state_n_plus()
                A = self.count_A(state_n, state_n_plus)
                return A * well_index[Components.WATER.value]
            else:
                return well_index[Components.OIL.value]
        else:
            return 0

    """
    def count_coefficient(self, flow, revert=False):
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()

        if not revert:
            main_cell = right_cell
            back_cell = left_cell
            pressure_main = left_cell.get_cell_state_n_plus().get_pressure_oil()
            pressure_back = right_cell.get_cell_state_n_plus().get_pressure_oil()
        else:
            main_cell = left_cell
            back_cell = right_cell
            pressure_main = right_cell.get_cell_state_n_plus().get_pressure_oil()
            pressure_back = left_cell.get_cell_state_n_plus().get_pressure_oil()

        if main_cell is None:
            return 0

        A = 1
        coeff = A * flow.get_water_flow() + flow.get_oil_flow()
        r_ost = -A * flow.get_water_flow() * (
            main_cell.get_cell_state_n_plus().get_pressure_cap() - back_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(main_cell.get_eq_index(), -r_ost - coeff * (pressure_main - pressure_back))
        return coeff
    """

    def count_b(self, flow):
        #========================== Стоит в матрице на (x, y + Nx * Ny). P(i+1, j, k) ==========================
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()
        if right_cell is None:
            return 0
        pressure_left = left_cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_right = right_cell.get_cell_state_n_plus().get_pressure_oil()

        state_n = right_cell.get_cell_state_n()
        state_n_plus = right_cell.get_cell_state_n_plus()
        A = self.count_A(state_n, state_n_plus)

        b = A * flow.get_water_flow() + flow.get_oil_flow()
        r_ost = -A * flow.get_water_flow() * (
            right_cell.get_cell_state_n_plus().get_pressure_cap() - left_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(right_cell.get_eq_index(), -r_ost - b * (pressure_left - pressure_right))
        return b

    def count_c(self, flow):
        #========================== Стоит в матрице на (x + Nx * Ny, y). P(i-1, j, k) ===============================
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()
        if left_cell is None:
            return 0
        pressure_left = left_cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_right = right_cell.get_cell_state_n_plus().get_pressure_oil()

        state_n = left_cell.get_cell_state_n()
        state_n_plus = left_cell.get_cell_state_n_plus()
        A = self.count_A(state_n, state_n_plus)
        c = A * flow.get_water_flow() + flow.get_oil_flow()
        r_ost = -A * flow.get_water_flow() * (
            left_cell.get_cell_state_n_plus().get_pressure_cap() - right_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(left_cell.get_eq_index(), -r_ost - c * (pressure_right - pressure_left))
        return c

    def count_d(self, flow):
        #========================== Стоит в матрице на (x + Nx * Ny, y). P(i, j + 1, k) ==========================
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()
        if right_cell is None:
            return 0
        pressure_left = left_cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_right = right_cell.get_cell_state_n_plus().get_pressure_oil()

        state_n = right_cell.get_cell_state_n()
        state_n_plus = right_cell.get_cell_state_n_plus()
        A = self.count_A(state_n, state_n_plus)
        d = A * flow.get_water_flow() + flow.get_oil_flow()
        r_ost = -A * flow.get_water_flow() * (
            right_cell.get_cell_state_n_plus().get_pressure_cap() - left_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(right_cell.get_eq_index(), -r_ost - d * (pressure_left - pressure_right))
        return d

    def count_e(self, flow):
        #========================== Стоит в матрице на (x + Nx * Ny, y). P(i, j - 1, k) ==========================
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()

        if left_cell is None:
            return 0

        state_n = left_cell.get_cell_state_n()
        state_n_plus = left_cell.get_cell_state_n_plus()
        A = self.count_A(state_n, state_n_plus)
        e = A * flow.get_water_flow() + flow.get_oil_flow()
        pressure_left = left_cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_right = right_cell.get_cell_state_n_plus().get_pressure_oil()
        r_ost = -A * flow.get_water_flow() * (
            left_cell.get_cell_state_n_plus().get_pressure_cap() - right_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(left_cell.get_eq_index(), -r_ost - e * (pressure_right - pressure_left))
        return e

    def count_g(self, flow):
        #========================== Стоит в матрице на (x + Nx * Ny, y). P(i, j, k + 1) ==========================
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()
        if right_cell is None:
            return 0
        pressure_left = left_cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_right = right_cell.get_cell_state_n_plus().get_pressure_oil()

        state_n = right_cell.get_cell_state_n()
        state_n_plus = right_cell.get_cell_state_n_plus()
        A = self.count_A(state_n, state_n_plus)
        g = A * flow.get_water_flow() + flow.get_oil_flow()
        if Layer.s_water_init != 0.0:
            gamma_water = (state_n_plus.get_s_water() * state_n_plus.get_ro_water() + left_cell.get_cell_state_n_plus().get_s_water() * left_cell.get_cell_state_n_plus().get_ro_water())\
                    / (state_n_plus.get_s_water() + left_cell.get_cell_state_n_plus().get_s_water())
        else:
            gamma_water = 0.0

        if Layer.s_oil_init != 0.0:
            gamma_oil = (
                      state_n_plus.get_s_oil() * state_n_plus.get_ro_oil() + left_cell.get_cell_state_n_plus().get_s_oil() * left_cell.get_cell_state_n_plus().get_ro_oil()) \
                      / (state_n_plus.get_s_oil() + left_cell.get_cell_state_n_plus().get_s_oil())
        else:
            gamma_oil = 0.0

        gravitation_water = (left_cell.z_coordinate - right_cell.z_coordinate) * gamma_water * Layer.g
        gravitation_oil = (left_cell.z_coordinate - right_cell.z_coordinate) * gamma_oil * Layer.g
        r_ost = -A * flow.get_water_flow() * (
            right_cell.get_cell_state_n_plus().get_pressure_cap() - left_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(right_cell.get_eq_index(), -r_ost - g * (pressure_left - pressure_right) + flow.get_oil_flow() * gravitation_oil + A * flow.get_water_flow() * gravitation_water)
        return g

    def count_f(self, flow):
        #========================== Стоит в матрице на (x + Nx * Ny, y). P(i, j, k - 1) ==========================
        left_cell = flow.get_left_cell()
        right_cell = flow.get_right_cell()
        if left_cell is None:
            return 0

        state_n = left_cell.get_cell_state_n()
        state_n_plus = left_cell.get_cell_state_n_plus()
        A = self.count_A(state_n, state_n_plus)
        f = A * flow.get_water_flow() + flow.get_oil_flow()
        if Layer.s_water_init != 0.0:
            gamma_water = (state_n_plus.get_s_water() * state_n_plus.get_ro_water() + right_cell.get_cell_state_n_plus().get_s_water() * right_cell.get_cell_state_n_plus().get_ro_water()) \
                / (state_n_plus.get_s_water() + right_cell.get_cell_state_n_plus().get_s_water())
        else:
            gamma_water = 0.0

        if Layer.s_oil_init != 0.0:
            gamma_oil = (state_n_plus.get_s_oil() * state_n_plus.get_ro_oil() + right_cell.get_cell_state_n_plus().get_s_oil() * right_cell.get_cell_state_n_plus().get_ro_oil()) \
                    / (state_n_plus.get_s_oil() + right_cell.get_cell_state_n_plus().get_s_oil())
        else:
            gamma_oil = 0.0

        gravitation_water = (right_cell.z_coordinate - left_cell.z_coordinate) * gamma_water * Layer.g
        gravitation_oil = (right_cell.z_coordinate - left_cell.z_coordinate) * gamma_oil * Layer.g
        pressure_left = left_cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_right = right_cell.get_cell_state_n_plus().get_pressure_oil()
        r_ost = -A * flow.get_water_flow() * (
            left_cell.get_cell_state_n_plus().get_pressure_cap() - right_cell.get_cell_state_n_plus().get_pressure_cap())
        self.solver_slau.add_nevyaz(left_cell.get_eq_index(), -r_ost - f * (pressure_right - pressure_left) + flow.get_oil_flow() * gravitation_oil + A * flow.get_water_flow() * gravitation_water)
        return f

    """
    def count_add_coefficients(self, cell_container, k, i, j):
        cell = cell_container.get_cell(k, i, j)
        pressure_n = cell.get_cell_state_n().get_pressure_oil()
        pressure_n_plus = cell.get_cell_state_n_plus().get_pressure_oil()
        pressure_cap_n = cell.get_cell_state_n().get_pressure_cap()

        self.b = np.append(self.b, self.count_b(
            cell.get_flow_coefficient('x', FlowComponents.plus.value)))  # Откусываем большой сдвиг с конца
        self.c = np.append(self.c, self.count_c(
            cell.get_flow_coefficient('x', FlowComponents.minus.value)))  # Откусываем большой сдвиг с начала
        self.d = np.append(self.d, self.count_d(
            cell.get_flow_coefficient('y', FlowComponents.plus.value)))  # Откусываем малый сдвиг с конца
        self.e = np.append(self.e, self.count_e(
            cell.get_flow_coefficient('y', FlowComponents.minus.value)))  # Откусываем малый сдвиг с начала
        self.g = np.append(self.g, self.count_g(
            cell.get_flow_coefficient('z', FlowComponents.plus.value)))  # Остается неизменным сдвиг + 1
        self.f = np.append(self.f, self.count_f(
            cell.get_flow_coefficient('z', FlowComponents.minus.value)))  # Остается неизменным сдвиг -1

        # Добавляем в середину
        well_oil = self.count_well_addition(cell)
        well_water = self.count_well_addition(cell, water=True)
        self.well_addition_oil = np.append(self.well_addition_oil, well_oil)  # Water False
        self.well_addition_water = np.append(self.well_addition_water, well_water)

        p_well = cell.well.p_well if cell.has_well else 0
        well_nevyaz = well_oil * (pressure_n_plus - p_well) + well_water * (pressure_n_plus - pressure_cap_n - p_well)
        self.coeff_well_nevyaz = np.append(self.coeff_well_nevyaz, well_nevyaz)

        a_d = self.count_a(cell)
        self.coeff_a_d = np.append(self.coeff_a_d, a_d)
        self.coeff_a_d_nevyaz = np.append(self.coeff_a_d_nevyaz, a_d * (pressure_n_plus - pressure_n))
        """

    """
    def end_coeff(self):
        self.solver_slau.add_vector_to_nevyaz(np.transpose([self.coeff_a_d_nevyaz]))
        self.solver_slau.add_vector_to_nevyaz(np.transpose([self.coeff_well_nevyaz]))

        self.a = - self.b - self.c - self.d - self.e - self.f - self.g - self.coeff_a_d - self.well_addition_oil - self.well_addition_water

        shift_Nx = self.solver_slau.get_shift_N_x()
        shift_Nxy = self.solver_slau.get_shift_N_xy()

        self.f = self.f[(shift_Nxy):]
        self.g = self.g[:-(shift_Nxy)]

        self.e = self.e[1:]
        self.d = self.d[:-1]

        self.c = self.c[(shift_Nx):]
        self.b = self.b[:(-shift_Nx)]

        diagonals = [self.a, self.e, self.d, self.c, self.b, self.f, self.g]
        shifts = [0, -1, 1, -shift_Nx, shift_Nx, -shift_Nxy, shift_Nxy]
        self.solver_slau.coefficient_matrix = diags(diagonals, shifts)
        # Для смотрения
        # smotrenie = self.solver_slau.coefficient_matrix.toarray()
        # print("be")
    """

    def generate_matrix(self, cell_container):
        matrix_size = Layer.N_x * Layer.N_y * Layer.N_z
        b = np.empty((0, matrix_size))  # Плюсовый
        c = np.empty((0, matrix_size))  # Минусовый
        d = np.empty((0, matrix_size))  # Плюсовый
        e = np.empty((0, matrix_size))  # Минусовый
        g = np.empty((0, matrix_size))  # Плюсовый
        f = np.empty((0, matrix_size))  # Минусовый
        well_addition_oil = np.empty((0, matrix_size))
        well_addition_water = np.empty((0, matrix_size))
        coeff_a_d = np.empty((0, matrix_size))
        coeff_a_d_nevyaz = np.empty((0, matrix_size))
        coeff_well_nevyaz = np.empty((0, matrix_size))

        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    cell = cell_container.get_cell(k, i, j)
                    pressure_n = cell.get_cell_state_n().get_pressure_oil()
                    pressure_n_plus = cell.get_cell_state_n_plus().get_pressure_oil()
                    pressure_cap_n = cell.get_cell_state_n().get_pressure_cap()

                    s_oil = cell.get_cell_state_n().get_s_oil()
                    s_water = cell.get_cell_state_n().get_s_water()

                    b = np.append(b, self.count_b(cell.get_flow_coefficient('x', FlowComponents.plus.value)))  # Откусываем большой сдвиг с конца
                    c = np.append(c, self.count_c(cell.get_flow_coefficient('x', FlowComponents.minus.value)))  # Откусываем большой сдвиг с начала
                    d = np.append(d, self.count_d(cell.get_flow_coefficient('y', FlowComponents.plus.value)))  # Откусываем малый сдвиг с конца
                    e = np.append(e, self.count_e(cell.get_flow_coefficient('y', FlowComponents.minus.value)))  # Откусываем малый сдвиг с начала
                    g = np.append(g, self.count_g(cell.get_flow_coefficient('z', FlowComponents.plus.value)))  # Остается неизменным сдвиг + 1
                    f = np.append(f, self.count_f(cell.get_flow_coefficient('z', FlowComponents.minus.value)))  # Остается неизменным сдвиг -1

                    # Добавляем в середину
                    well_oil = self.count_well_addition(cell)
                    well_water = self.count_well_addition(cell, water=True)
                    well_addition_oil = np.append(well_addition_oil, well_oil)  # Water False
                    well_addition_water = np.append(well_addition_water, well_water)

                    # TODO: спросить какое давление капплярное сюда пихать
                    p_well = cell.well.p_well if cell.has_well else 0

                    p_average = pressure_n_plus - s_water*pressure_cap_n
                    well_nevyaz = well_oil * (p_average - p_well) + well_water * (p_average - p_well)
                    coeff_well_nevyaz = np.append(coeff_well_nevyaz, well_nevyaz)

                    a_d = self.count_a(cell)
                    coeff_a_d = np.append(coeff_a_d, a_d)
                    coeff_a_d_nevyaz = np.append(coeff_a_d_nevyaz, a_d * (pressure_n_plus - pressure_n))

        # Невязка a_d
        self.solver_slau.add_vector_to_nevyaz(np.transpose([coeff_a_d_nevyaz]))
        self.solver_slau.add_vector_to_nevyaz(np.transpose([coeff_well_nevyaz]))

        a = - b - c - d - e - f - g - coeff_a_d - well_addition_oil - well_addition_water

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
        #smotrenie = self.solver_slau.coefficient_matrix.toarray()
        #print("be")

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
                    state_n_plus.set_pressure_water(state_n_plus.get_pressure_oil() - state_n_plus.get_pressure_cap())
                    eq_index += 1

    def update_saturation(self, cell_container):

        def count_coeff_2(flow_plus, flow_minus):
            right_state_plus = flow_plus.get_right_cell().get_cell_state_n_plus()
            left_state_plus = flow_plus.get_left_cell().get_cell_state_n_plus()

            right_state_minus = flow_minus.get_right_cell().get_cell_state_n_plus()
            left_state_minus = flow_minus.get_left_cell().get_cell_state_n_plus()

            coeff_2_part = flow_plus.get_water_flow() * (right_state_plus.get_pressure_oil() - left_state_plus.get_pressure_oil())\
                           + flow_plus.get_water_flow() * (right_state_plus.get_pressure_cap() - left_state_plus.get_pressure_cap())\
                           + flow_minus.get_water_flow() * (left_state_minus.get_pressure_oil() - right_state_minus.get_pressure_oil())\
                           + flow_minus.get_water_flow() * (left_state_minus.get_pressure_cap() - right_state_minus.get_pressure_cap())

            # TODO: заменить расчет плотности в gamma если нужно
            gamma = left_state_plus.get_ro_water() * Layer.g
            gamma_part = flow_plus.get_water_flow() * gamma * (flow_plus.get_right_cell().z_coordinate - flow_plus.get_left_cell().z_coordinate)\
                         + flow_minus.get_water_flow() * gamma * (flow_minus.get_left_cell().z_coordinate - flow_minus.get_right_cell().z_coordinate)

            coeff_2_part -= gamma_part
            return coeff_2_part

        for k in range(1, Layer.N_z - 1):
            for i in range(1, Layer.N_x - 1):
                for j in range(1, Layer.N_y - 1):
                    #==================== Все значения для воды ========================
                    cell = cell_container.get_cell(k, i, j)
                    flow_x_plus = cell.get_flow_coefficient('x', FlowComponents.plus.value)
                    flow_x_minus = cell.get_flow_coefficient('x', FlowComponents.minus.value)

                    flow_y_plus = cell.get_flow_coefficient('y', FlowComponents.plus.value)
                    flow_y_minus = cell.get_flow_coefficient('y', FlowComponents.minus.value)

                    flow_z_plus = cell.get_flow_coefficient('z', FlowComponents.plus.value)
                    flow_z_minus = cell.get_flow_coefficient('z', FlowComponents.minus.value)

                    state_n = cell.get_cell_state_n()
                    state_n_plus = cell.get_cell_state_n_plus()
                    #================ water ================
                    # TODO: Запихать эти штуки класс для воды или нефти
                    ro_der = Layer.water.ro_water_0 * Layer.water.c_f_water
                    fi_der = Layer.fi_0 * Layer.c_r

                    p_cap_der = Layer.count_p_cap_graph_der(state_n.get_s_water())
                    pressure_water = state_n_plus.get_pressure_water()

                    #================ коэффициенты d ===============
                    d_11 = (1.0 / self.tau) * state_n.get_s_water() * (
                    state_n.get_fi() * ro_der + fi_der * state_n_plus.get_ro_water())
                    d_12 = (1.0 / self.tau) * (
                    state_n_plus.get_fi() * state_n_plus.get_ro_water() - state_n.get_s_water() * state_n.get_fi() * state_n_plus.get_ro_water() * p_cap_der)

                    coeff_2 = count_coeff_2(flow_x_plus, flow_x_minus) + count_coeff_2(flow_y_plus, flow_y_minus) + count_coeff_2(flow_z_plus, flow_z_minus)

                    coeff_3 = -(d_11 / d_12) * (state_n_plus.get_pressure_oil() - state_n.get_pressure_oil())

                    s_water_new = state_n.get_s_water() + (1.0 / d_12) * coeff_2 + coeff_3
                    if cell.has_well:
                        s_water_new = s_water_new - (1.0 / d_12) * cell.well.get_water_well_index() * (pressure_water - cell.well.p_well)
                    state_n_plus.set_s_water(s_water_new)
                    state_n_plus.set_s_oil(1.0 - s_water_new)

    def update_pressure_cap(self, cell_container):
        for k in range(1, Layer.N_z - 1):
            for i in range(1, Layer.N_x - 1):
                for j in range(1, Layer.N_y - 1):
                    cell = cell_container.get_cell(k, i, j)
                    s_water = cell.get_cell_state_n_plus().get_s_water()
                    p_cap = Layer.count_pressure_cap(s_water)
                    cell.get_cell_state_n_plus().set_pressure_cap(p_cap)

    def check_pressure_convergence(self, cell_container, delta_k):
        conv = []
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    cell = cell_container.get_cell(k, i, j)
                    pressure_oil = cell.get_cell_state_n_plus().get_pressure_oil()
                    conv.append(delta_k[cell.get_eq_index()] / pressure_oil)
        norm = self.count_norm(delta_k)
        if self.count_norm(delta_k) > 0.1:
            print(f'Pressure norm {norm} > 10%')
            return True
        else:
            return False

    def count_debit(self, cell_container):
        for well_index in Layer.wells:
            cell = cell_container.get_cell(well_index[0], well_index[1], well_index[2])
            pressure = cell.get_cell_state_n_plus().get_pressure_water()

    def main_cycle(self, cell_container):
        #self.set_zero_coeffs()
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    self.recount_property(cell_container, k, i, j)
                    self.count_cell_flow(cell_container, k, i, j)
                    #self.count_add_coefficients(cell_container, k, i, j)
        #self.end_coeff()
