from .Layer import Layer


class ThreeDimOilWaterImpes:
    def __init__(self):
        self.tau_default = 86400
        self.tau = self.tau_default

    def generate_delta_k(self):
        # TODO: Сгенерировать начальное приближение
        pass

    def count_norm(self, delta_list_k):
        # TODO: Посчитать нормуд дельта как модуль всех элементов ты поняла
        pass

    def check_norm(self, delta_list_k):
        # TODO: Сравнить ее с максимальным значением. Пропихнуть его при инициализации кстати да
        pass

    def recount_properties(self, cell_container):
        # TODO: Пересчитаь все параметры типо давления и прочего
        pass

    def count_flows(self, flows):
        # TODO: Проинициализировать значения всех потоков
        pass

    def generate_matrix(self):
        # TODO: Собственно сгенерировать всю матрицу. Ну которая семидиагональная ага да
        pass

    def solve_slau(self):
        # TODO: Реши систему. Тащемта можешь использовать стандартный решатель
        pass

    def update_pressure(self, cell_container, delta_k):
        # TODO: Обнови давление ага да
        pass

    def update_saturation(self, cell_container, flows):
        # TODO: здесь обнови насыщенность
        pass

    def update_pressure_cap(self, cell_container):
        # TODO: КАПИЛЯРКА АГА (обнови) Перестань общаться с собой через комментарии к коду!!
        pass

    def check_pressure_convergence(self, cell_container, delta_k):
        # TODO: Проверить сходимость по давлению
        pass