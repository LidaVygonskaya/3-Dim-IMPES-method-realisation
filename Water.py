from .Layer import Layer


class Water:
    def __init__(self):
        # TODO: изменить значение P_01
        self.ro_water_0 = 1000.0
        self.c_f_water = (10.0 ** (-4)) / Layer.atm
        self.P_01 = Layer.P_01

    def count_ro(self, pressure_water):
        return self.ro_water_0 * (1.0 + self.c_f_water * (pressure_water - self.P_01))

    @staticmethod
    def count_k_r(s_water):
        return s_water ** 2.0