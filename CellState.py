import numpy as np

from Flow import Flow
from Layer import Layer
from Enums import Components


class CellState:

    def __init__(self):
        self.pressure_oil_water = np.zeros(Layer.components_count)  # [P_oil, P_water]
        self.pressure_cap = 0.0
        self.s_oil_water = np.zeros(Layer.components_count)  # [S_oil, S_water]
        self.ro_oil_water = np.zeros(Layer.components_count)  # [ro_oil, ro_water]
        self.k_r_oil_water = np.zeros(Layer.components_count)  # [k_r_oil, k_r_water]

        #self.t_minus_oil_water = np.zeros(Layer.components_count)  # [t_minus_oil, t_minus_water]
        #self.t_plus_oil_water = np.zeros(Layer.components_count)  # [t_plus_oil, t_plus_water]

        self.fi = 0.0

        # TODO: Может перенести расчет этого говна в другое место
        self.c1_p = 0.0
        self.c2_p = 0.0



    #Getters
    def get_pressure_oil(self):
        return self.pressure_oil_water[Components.OIL.value]  # Вернет конкретно значение

    def get_pressure_water(self):
        return self.pressure_oil_water[Components.WATER.value]

    def get_pressure_cap(self):
        return self.pressure_cap

    def get_components_pressure(self):
        return self.pressure_oil_water

    def get_components_saturation(self):
        return self.s_oil_water

    def get_s_water(self):
        return self.s_oil_water[Components.WATER.value]

    def get_s_oil(self):
        return self.s_oil_water[Components.OIL.value]

    def get_ro_water(self):
        return self.ro_oil_water[Components.WATER.value]

    def get_ro_oil(self):
        return self.ro_oil_water[Components.OIL.value]

    def get_components_ro(self):
        return self.ro_oil_water

    def get_components_k_r(self):
        return self.k_r_oil_water

    def get_k_r_oil(self):
        return self.k_r_oil_water[Components.OIL.value]

    def get_k_r_water(self):
        return self.k_r_oil_water[Components.WATER.value]

    def get_fi(self):
        return self.fi

    def get_c1_p(self):
        pass

    def get_c2_p(self):
        pass


    #Setters
    # Здесь устанавливаем значение тупо по значениям
    def set_pressure_oil(self, pressure_oil):
        self.pressure_oil_water[Components.OIL.value] = pressure_oil

    def set_pressure_water(self, pressure_water):
        self.pressure_oil_water[Components.WATER.value] = pressure_water

    def set_s_water(self, s_water):
        self.s_oil_water[Components.WATER.value] = s_water

    def set_s_oil(self, s_oil):
        self.s_oil_water[Components.OIL.value] = s_oil

    def set_ro_water(self, ro_water):
        self.ro_oil_water[Components.WATER.value] = ro_water

    def set_ro_oil(self, ro_oil):
        self.ro_oil_water[Components.OIL.value] = ro_oil

    def set_k_r(self, k_r, index):
        self.k_r_oil_water[index] = k_r

    def set_ro(self, ro, index):
        self.ro_oil_water[index] = ro

    def set_fi(self, fi):
        self.fi = fi

    def set_pressure_cap(self, pressure_cap):
        self.pressure_cap = pressure_cap


    # TODO : Перенести данные коэффициенты
    def set_c1_p(self, c1_p):
        self.c1_p = c1_p

    def set_c2_p(self, c2_p):
        self.c2_p = c2_p


    #Set n to n plus 1
    def set_equals_to(self, cell_state_new):
        self.pressure_oil_water[Components.OIL.value] = cell_state_new.pressure_oil_water[Components.OIL.value]
        self.pressure_oil_water[Components.WATER.value] = cell_state_new.pressure_oil_water[Components.WATER.value]
        self.pressure_cap = cell_state_new.pressure_cap
        self.s_oil_water[Components.OIL.value] = cell_state_new.s_oil_water[Components.OIL.value]
        self.s_oil_water[Components.WATER.value] = cell_state_new.s_oil_water[Components.WATER.value]
        self.ro_oil_water[Components.OIL.value] = cell_state_new.ro_oil_water[Components.OIL.value]
        self.ro_oil_water[Components.WATER.value] = cell_state_new.ro_oil_water[Components.WATER.value]
        self.k_r_oil_water[Components.OIL.value] = cell_state_new.k_r_oil_water[Components.OIL.value]
        self.k_r_oil_water[Components.WATER.value] = cell_state_new.k_r_oil_water[Components.WATER.value]
        self.fi = cell_state_new.fi

        # TODO: Перенести данные вещи в построение матрицы
        self.c1_p = cell_state_new.c1_p
        self.c2_p = cell_state_new.c2_p