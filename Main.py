import numpy as np

from .ThreeDimOilWaterImpes import ThreeDimOilWaterImpes
from .CellContainer import CellContainer
from .Flow import Flow

impes = ThreeDimOilWaterImpes()
# TODO: запихать сюда еще и решателя (а может не надо???)

cell_container = CellContainer()  # Проверь на счет eq_index. Внутри реализации написано чо каво
cell_container.initialize_cells()

flows = Flow.initialize_flow_array()

time = impes.tau  # Сразу обозначим это как первый шаг по времени, потому что нулевой у нас есть


# TODO: А вот здесь добавь шаг по времени. Пока посчитаем только для одного
delta_k = impes.generate_delta_k()  # Генерируем начальное приближение (Я больше не лист я куб нах)
cell_container.equate_cell_states()  # State_n = State_n_plus

while impes.check_norm(delta_k):
    # solver.set_zero()
    impes.recount_properties(cell_container)
    impes.count_flows(flows)
    impes.generate_matrix()

    impes.solve_slau()
    #delta_k = solver_slau.get_result() Не знаю, Лидос, нужны ли тебе эти строчки, смотри сама
    #solver_slau.clear_result()

    if impes.check_pressure_convergence(cell_container, delta_k):
        impes.tau = impes.tau / 2.0
        cell_container.equate_cell_states(to_previous=True)
    else:
        impes.update_pressure(cell_container, delta_k)

impes.recount_properties(cell_container)
impes.count_flows(flows)
impes.update_saturation(cell_container, flows)
impes.update_pressure_cap(cell_container)