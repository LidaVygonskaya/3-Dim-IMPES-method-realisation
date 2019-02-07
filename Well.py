import math
import numpy as np

from Layer import Layer
from Oil import Oil
from Water import Water
from Enums import Components


class Well:
    def __init__(self):
        extractive = self.define_well_type()
        self.p_well = Layer.P_well_extractive if extractive else Layer.P_well_delivery
        self.r_well = Layer.r_well

        # TODO: менять пересчет насыщенности в зависимости от типа скважины
        self.s_well_water = Layer.s_well_water
        self.s_well_oil = Layer.s_well_oil

        self.well_index_oil_water = np.zeros(Layer.components_count)
        self.count_well_index()

    def count_well_index(self):
        k = self.count_k()
        re = self.count_re()
        k_r_oil = Oil.count_k_r(self.s_well_oil)
        k_r_water = Water.count_k_r(self.s_well_water)
        well_oil = 2.0 * math.pi * k * Layer.h_z * k_r_oil / (math.log(re / self.r_well) * Layer.mu_oil)
        well_water = 2.0 * math.pi * k * Layer.h_z * k_r_water / (math.log(re / self.r_well) * Layer.mu_water)
        self.well_index_oil_water[Components.OIL.value] = well_oil / Layer.V_ijk
        self.well_index_oil_water[Components.WATER.value] = well_water / Layer.V_ijk

    def count_k(self):
        return math.sqrt(Layer.k_x * Layer.k_y)

    def count_re(self):
        nominator = math.sqrt(Layer.k_y / Layer.k_x) * (Layer.h_x) ** 2 + math.sqrt(Layer.k_x / Layer.k_y) * (Layer.h_y) ** 2
        denominator = math.pow(Layer.k_y / Layer.k_x, 0.25) + math.pow(Layer.k_x / Layer.k_y, 0.25)
        return 0.28 * math.sqrt(nominator) / denominator

    def define_well_type(self):
        # TODO: дописать эту функцию
        extractive = False
        p_average = Layer.pressure_water_init * Layer.s_water_init + Layer.pressure_oil_init * Layer.s_oil_init
        return extractive