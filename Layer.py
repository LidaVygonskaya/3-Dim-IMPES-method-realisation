import numpy as np

from Oil import Oil
from Water import Water


class Layer:
    atm = 101325.0
    p_cap_filename = 'Pcap(Sw).txt'

    components_count = 2  # Количество компонент
    dim = 3  # Размерность пространства

    # Количество ячеек
    N_x = 7
    N_y = 7
    N_z = 7

    # Координаты и шаги
    # TODO: По разной координатной оси разный шаг по пространству
    # TODO: Сейчас пока что у тебя один шаг. Надо будет все поменять
    x_0 = 0.0
    x_N = 100.0
    h = (x_N - x_0) / (N_x - 1)
    h_x = (x_N - x_0) / (N_x - 1)
    h_y = (x_N - x_0) / (N_y - 1)
    h_z = (x_N - x_0) / (N_z - 1)
    h_array = [h_x, h_y, h_z]

    s_water_init = 0.0  # Начальная насыщенность воды
    s_oil_init = 1.0  # Начальная насыщенность нефти. Вычисляем через воду конечно, но пусть будет на всякий

    pressure_oil_init = 80 * atm  # Начальное давление нефти
    pressure_water_init = 80 * atm  # Начальное давление воды. Вычисляется через капилярку, но пусть будет на всякий

    fi_0 = 0.2
    c_r = (10.0 ** (-5)) / atm

    P_01 = Water.P_01
    P_02 = Oil.P_02

    k = (9.868233 * (10 ** (-13))) * 10 ** (-3)
    # TODO: Создать коэффициенты проницаемости
    k_x = 0.0
    k_y = k
    k_z = 0.0
    k_array = [k_x, k_y, k_z]
    mu_oil = 10.0 * (10.0 ** (-3))
    mu_water = 10.0 ** (-3)
    mu_oil_water = [mu_oil, mu_water]
    oil = Oil()
    water = Water()
    components = [oil, water]

    delta_0 = 1000.0  # Начальная прибавка по давлению в паскалях

    @staticmethod
    def get_h(direction):
        if direction == 'x':
            return Layer.h_array[0]
        if direction == 'y':
            return Layer.h_array[1]
        if direction == 'z':
            return Layer.h_array[2]

    @staticmethod
    def count_pressure_cap(s_water):
        # TODO: Организовать подачу вектора на вход. Нельзя забывать что теперь это уже не скаляр
        # TODO: Нужно пробегаться по каждому значению x y z насыщенности в массиве
        pressure_cap_graph = {}  # s_water: pressure_cap
        file = open(Layer.p_cap_filename, 'r')
        for line in file.readlines():
            line = line.rstrip()
            line = line.split('\t')
            pressure_cap_graph.update({float(line[0]): float(line[1]) * Layer.atm})
        s_w_graph = list(pressure_cap_graph.keys())
        for i in range(len(s_w_graph)):
            if s_water <= s_w_graph[i]:
                p_cap = pressure_cap_graph.get(s_w_graph[i - 1]) + (pressure_cap_graph.get(s_w_graph[i])
                                                                    - pressure_cap_graph.get(s_w_graph[i - 1])) \
                                                                    / (s_w_graph[i] - s_w_graph[i - 1]) \
                                                                    * (s_water - s_w_graph[i - 1])
                # p_cap = 0
                break
        return p_cap

    @classmethod
    def count_fi(cls,pressure_oil):
        return cls.fi_0 * (1.0 + cls.c_r * (pressure_oil - cls.P_01))
