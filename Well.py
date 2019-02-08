import math
import numpy as np

from Layer import Layer
from Oil import Oil
from Water import Water
from Enums import Components


class Well:
    def __init__(self, cell):
        self.p_well = Layer.P_well_extractive if Layer.productive else Layer.P_well_delivery
        self.r_well = Layer.r_well

        self.s_well_water = cell.get_cell_state_n_plus().get_s_water() if Layer.productive else Layer.s_well_water
        self.s_well_oil = cell.get_cell_state_n_plus().get_s_oil() if Layer.productive else Layer.s_well_oil

        self.well_index_oil_water = np.zeros(Layer.components_count)
        self.count_well_index(cell)

    def count_well_index(self, cell, horizontal=False):
        k = self.count_k(horizontal)
        re = self.count_re(horizontal)
        k_r_oil = Oil.count_k_r(self.s_well_oil)
        k_r_water = Water.count_k_r(self.s_well_water)
        delta_dim = Layer.h_y if horizontal else Layer.h_z
        well_oil = 2.0 * math.pi * k * delta_dim * k_r_oil / (math.log(re / self.r_well) * Layer.mu_oil)
        well_water = 2.0 * math.pi * k * delta_dim * k_r_water / (math.log(re / self.r_well) * Layer.mu_water)
        ro_oil = Oil.count_ro(cell.get_cell_state_n_plus().get_pressure_water()) if Layer.productive else Oil.count_ro(self.p_well)
        ro_water = cell.get_cell_state_n_plus().get_ro_water()
        self.well_index_oil_water[Components.OIL.value] = well_oil * ro_oil / Layer.V_ijk
        self.well_index_oil_water[Components.WATER.value] = well_water * ro_water / Layer.V_ijk

    def get_water_well_index(self):
        return self.well_index_oil_water[Components.WATER.value]

    def get_oil_well_index(self):
        return self.well_index_oil_water[Components.OIL.value]

    def count_k(self, horizontal=False):
        if horizontal:
            return math.sqrt(Layer.k_x * Layer.k_z)
        return math.sqrt(Layer.k_x * Layer.k_y)

    def count_re(self, horizontal=False):
        if horizontal:
            nominator = math.sqrt(Layer.k_y / Layer.k_x) * (Layer.h_x) ** 2 + math.sqrt(Layer.k_x / Layer.k_y) * (Layer.h_y) ** 2
            denominator = math.pow(Layer.k_y / Layer.k_x, 0.25) + math.pow(Layer.k_x / Layer.k_y, 0.25)
        else:
            nominator = math.sqrt(Layer.k_y / Layer.k_x) * (Layer.h_x) ** 2 + math.sqrt(Layer.k_x / Layer.k_y) * (Layer.h_y) ** 2
            denominator = math.pow(Layer.k_y / Layer.k_x, 0.25) + math.pow(Layer.k_x / Layer.k_y, 0.25)
        return 0.28 * math.sqrt(nominator) / denominator

    def define_well_type(self):
        # TODO: дописать эту функцию
        extractive = False
        p_average = Layer.pressure_water_init * Layer.s_water_init + Layer.pressure_oil_init * Layer.s_oil_init
        return extractive
