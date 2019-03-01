import math
import numpy as np

from Layer import Layer
from Oil import Oil
from Water import Water
from Enums import Components, Boundary


class Well:
    def __init__(self, cell, well_index, horizontal):
        self.p_well = Layer.P_well_extractive if Layer.productive else Layer.P_well_delivery
        self.r_well = Layer.r_well

        self.s_well_water = cell.get_cell_state_n_plus().get_s_water() if Layer.productive else Layer.s_well_water
        self.s_well_oil = cell.get_cell_state_n_plus().get_s_oil() if Layer.productive else Layer.s_well_oil

        self.well_index_oil_water = np.zeros(Layer.components_count)
        self.count_well_index(cell, well_index, horizontal)

    def count_well_index(self, cell, well_index,  horizontal=False):
        k = self.count_k(horizontal)
        re = self.count_re(well_index, horizontal=horizontal)
        k_r_oil = Oil.count_k_r(self.s_well_oil)
        k_r_water = Water.count_k_r(self.s_well_water)
        delta_dim = Layer.h_y if horizontal else Layer.h_z
        well_oil = 2.0 * math.pi * k * delta_dim * k_r_oil / (math.log(re / self.r_well) * Layer.mu_oil)
        well_water = 2.0 * math.pi * k * delta_dim * k_r_water / (math.log(re / self.r_well) * Layer.mu_water)
        ro_oil = Oil.count_ro(cell.get_cell_state_n_plus().get_pressure_oil()) if Layer.productive else Oil.count_ro(self.p_well)
        ro_water = Water.count_ro(cell.get_cell_state_n_plus().get_pressure_oil()) if Layer.productive else Oil.count_ro(self.p_well)
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

    def count_re(self, well_index, horizontal=False):
        if horizontal:
            nominator = math.sqrt(Layer.k_y / Layer.k_z) * (Layer.h_z) ** 2 + math.sqrt(Layer.k_z / Layer.k_y) * (Layer.h_y) ** 2
            denominator = math.pow(Layer.k_y / Layer.k_z, 0.25) + math.pow(Layer.k_z / Layer.k_y, 0.25)
        else:
            nominator = math.sqrt(Layer.k_y / Layer.k_x) * (Layer.h_x) ** 2 + math.sqrt(Layer.k_x / Layer.k_y) * (Layer.h_y) ** 2
            denominator = math.pow(Layer.k_y / Layer.k_x, 0.25) + math.pow(Layer.k_x / Layer.k_y, 0.25)
        coeff = self.define_coefficient(well_index)
        return coeff * math.sqrt(nominator) / denominator

    def define_well_type(self):
        # TODO: дописать эту функцию
        extractive = False
        p_average = Layer.pressure_water_init * Layer.s_water_init + Layer.pressure_oil_init * Layer.s_oil_init
        return extractive

    def boundary(self, well_index):
        k, i, j = well_index
        count = 0

        if i == 0:
            count += 1
        elif i == Layer.N_x - 1:
            count += 1

        if j == 0:
            count += 1
        elif j == Layer.N_y - 1:
            count += 1

        if k == 0:
            count += 1
        elif k == Layer.N_z - 1:
            count += 1
        return count

    def define_coefficient(self, well_index):
        count = self.boundary(well_index)

        f_geo = 0.24879
        w_frac = 1
        if count == Boundary.BOUNDARY.value:
            f_geo = 0.235
            w_frac = 1
        elif count == Boundary.EDGE.value:
            f_geo = 0.229
            w_frac = 1

        return 2.0 * f_geo / (math.sqrt(math.pi * w_frac))
