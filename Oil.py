from .Layer import Layer


class Oil:
    def __init__(self):
        # TODO: изменить значение P_02
        self.ro_oil_0 = 900.0
        self.c_f_oil = (10.0 ** (-4)) / Layer.atm
        self.P_02 = Layer.P_02

    def count_ro(self, pressure_oil):
        return self.ro_oil_0 * (1.0 + self.c_f_oil * (pressure_oil - self.P_02))

    @staticmethod
    def count_k_r(s_oil):
        return s_oil ** 2.0