class Water:
    atm = 101325.0
    P_01 = 80 * atm
    ro_water_0 = 1000.0
    c_f_water = (10.0 ** (-4)) / atm

    def __init__(self):
        # TODO: изменить значение P_01
        self.ro_water_0 = 1000.0
        self.c_f_water = (10.0 ** (-4)) / Water.atm

    @staticmethod
    def count_ro(pressure_water):
        return Water.ro_water_0 * (1.0 + Water.c_f_water * (pressure_water - Water.P_01))

    @staticmethod
    def count_k_r(s_water):
        return s_water ** 2.0