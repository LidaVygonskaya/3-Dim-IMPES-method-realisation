class Oil:
    atm = 101325.0
    P_02 = 80 * atm

    def __init__(self):
        # TODO: изменить значение P_02
        self.ro_oil_0 = 900.0
        self.c_f_oil = 0 #(10.0 ** (-4)) / Oil.atm

    def count_ro(self, pressure_oil):
        return self.ro_oil_0 * (1.0 + self.c_f_oil * (pressure_oil - Oil.P_02))

    @staticmethod
    def count_k_r(s_oil):
        return s_oil ** 2.0