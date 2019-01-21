import numpy as np

from ThreeDimOilWaterImpes import ThreeDimOilWaterImpes
from CellContainer import CellContainer
from Flow import Flow
from SolverSlau import SolverSlau

solver_slau = SolverSlau()
impes = ThreeDimOilWaterImpes(solver_slau)

cell_container = CellContainer()  # Проверь на счет eq_index. Внутри реализации написано чо каво
cell_container.initialize_cells()
cell = cell_container.get_cell(3, 3, 3)

cell_container.get_cell(3, 3, 3).get_cell_state_n().set_pressure_oil(100*101325)
cell_container.get_cell(3, 3, 3).get_cell_state_n_plus().set_pressure_oil(100*101325 + 1000)

cell_container.initialize_flows()
Flow.initialize_flow(cell_container)
#flows = Flow.initialize_flow_array(cell_container)

time = impes.tau  # Сразу обозначим это как первый шаг по времени, потому что нулевой у нас есть

# TODO: А вот здесь добавь шаг по времени. Пока посчитаем только для одного
delta_k = impes.generate_delta_k()  # Генерируем начальное приближение (Я больше не лист я куб нах) Не нихуя все таки тупой массив
#cell_container.equate_cell_states()  # State_n = State_n_plus

while impes.check_norm(delta_k):
    impes.recount_properties(cell_container)
    impes.count_cells_flows(cell_container)
    impes.generate_matrix(cell_container)

    impes.solve_slau()
    delta_k = impes.solver_slau.get_result()
    f = open('govno2.txt', 'w')
    for de in delta_k:
        f.write(str(de) + '\n')
    f.close()
    impes.solver_slau.clear_result()
    impes.update_pressure(cell_container, delta_k)

    # TODO: проверять сходимость если нужно опять же я хз
    #delta_k = solver_slau.get_result() Не знаю, Лидос, нужны ли тебе эти строчки, смотри сама
    #solver_slau.clear_result()

    #if impes.check_pressure_convergence(cell_container, delta_k):
    #    impes.tau = impes.tau / 2.0
    #    cell_container.equate_cell_states(to_previous=True)
    #else:
    #    impes.update_pressure(cell_container, delta_k)

# TODO: Построить распределение давлений


impes.recount_properties(cell_container)
impes.count_flows(flows)
impes.update_saturation(cell_container, flows)
impes.update_pressure_cap(cell_container)