import numpy as np


class Layer:
    atm = 101325.0
    p_cap_filename = 'govno.txt'

    components_count = 2  # Количество компонент
    dim = 3  # Размерность пространства

    # Количество ячеек
    N_x = 100
    N_y = 100
    N_z = 100

    s_water_init = 0.5  # Начальная насыщенность воды
    s_oil_init = 0.5  # Начальная насыщенность нефти. Вычисляем через воду конечно, но пусть будет на всякий

    pressure_oil_init = 200 * atm  # Начальное давление нефти
    pressure_water_init = 200 * atm  # Начальное давление воды. Вычисляется через капилярку, но пусть будет на всякий

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
