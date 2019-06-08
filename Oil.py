class Oil:
    atm = 101325.0
    P_02 = 80 * atm
    ro_oil_0 = 900.0
    c_f_oil = (10.0 ** (-4)) / atm

    def __init__(self):
        self.ro_oil_0 = 900.0
        self.c_f_oil = (10.0 ** (-4)) / Oil.atm

    @staticmethod
    def count_ro(pressure_oil):
        return Oil.ro_oil_0 * (1.0 + Oil.c_f_oil * (pressure_oil - Oil.P_02))

    @staticmethod
    def count_k_r(s_oil):
        return s_oil ** 2.0