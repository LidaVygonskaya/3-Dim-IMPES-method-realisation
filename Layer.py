from Oil import Oil
from Water import Water


class Layer:
    atm = 101325.0
    g = 9.80665
    p_cap_filename = 'Pcap(Sw)_linear.txt'
    folder = '2D_A_2'

    components_count = 2  # Количество компонент
    dim = 3  # Размерность пространства

    # Количество ячеек
    N_x = 21
    N_y = 21
    N_z = 1

    # Координаты и шаги
    # TODO: По разной координатной оси разный шаг по пространству
    # TODO: Сейчас пока что у тебя один шаг. Надо будет все поменять
    x_0 = 0.0
    x_N = 200.0
    y_N = 200.0
    z_N = 1.0
    h = (x_N - x_0) / (N_x - 1)
    h_x = (x_N - x_0) / (N_x - 1)
    h_y = (y_N - x_0) / (N_y - 1)
    h_z = (z_N - x_0) #/ (N_z - 1)
    h_array = [h_x, h_y, h_z]
    V_ijk = h_x * h_y * h_z

    s_water_init = 0.5  # Начальная насыщенность воды
    s_oil_init = 1.0 - s_water_init  # Начальная насыщенность нефти. Вычисляем через воду конечно, но пусть будет на всякий

    pressure_oil_init = 80 * atm  # Начальное давление нефти
    pressure_water_init = 80 * atm  # Начальное давление воды. Вычисляется через капилярку, но пусть будет на всякий

    fi_0 = 0.2
    c_r = (10.0 ** (-5)) / atm

    P_01 = Water.P_01
    P_02 = Oil.P_02

    k = (9.868233 * (10 ** (-13))) * 10 ** (-2)
    # TODO: Создать коэффициенты проницаемости
    k_x = k
    k_y = k
    k_z = k
    k_array = [k_x, k_y, k_z]
    mu_oil = 10.0 * (10.0 ** (-3))
    mu_water = 10.0 * (10.0 ** (-3))
    mu_oil_water = [mu_oil, mu_water]
    oil = Oil()
    water = Water()
    components = [oil, water]

    delta_0 = 1000.0  # Начальная прибавка по давлению в паскалях

    # ========================== Скважинка ===========================
    P_well_delivery = 130 * atm  # injection
    P_well_extractive = 30 * atm  # production
    r_well = 0.108
    s_well_water = 1.0
    s_well_oil = 1.0 - s_well_water

    # k, i, j
    #wells = [(0, k, 100) for k in range(45, 56)]
    wells = [(0, 11, 11)]

    productive = True
    horizontal = False
    # ================================================================
    @staticmethod
    def get_h(direction):
        if direction == 'x':
            return Layer.h_array[0]
        if direction == 'y':
            return Layer.h_array[1]
        if direction == 'z':
            return Layer.h_array[2]

    @staticmethod
    def get_pressure_cap_graph():
        pressure_cap_graph = {}  # s_water: pressure_cap
        file = open(Layer.p_cap_filename, 'r')
        for line in file.readlines():
            line = line.rstrip()
            line = line.split('\t')
            pressure_cap_graph.update({float(line[0]): float(line[1]) * Layer.atm})
        return pressure_cap_graph


    @staticmethod
    def count_pressure_cap(s_water):
        pressure_cap_graph = Layer.get_pressure_cap_graph()  # s_water: pressure_cap
        s_w_graph = list(pressure_cap_graph.keys())
        for i in range(len(s_w_graph)):
            if s_water <= s_w_graph[i]:
                p_cap = pressure_cap_graph.get(s_w_graph[i - 1]) + (pressure_cap_graph.get(s_w_graph[i])
                                                                    - pressure_cap_graph.get(s_w_graph[i - 1])) \
                                                                    / (s_w_graph[i] - s_w_graph[i - 1]) \
                                                                    * (s_water - s_w_graph[i - 1])
                break
        return p_cap

    @staticmethod
    def count_p_cap_graph_der(s_water):
        pressure_cap_graph = Layer.get_pressure_cap_graph()
        s_w_graph = list(pressure_cap_graph.keys())
        for i in range(len(s_w_graph)):
            if s_water <= s_w_graph[i]:
                p_cap_der = (pressure_cap_graph.get(s_w_graph[i]) - pressure_cap_graph.get(s_w_graph[i - 1])) / (
                s_w_graph[i] - s_w_graph[i - 1])
                break
        return p_cap_der

    @classmethod
    def count_fi(cls, pressure_oil):
        return cls.fi_0 * (1.0 + cls.c_r * (pressure_oil - cls.P_01))