import numpy as np

from Layer import Layer
from ThreeDimOilWaterImpes import ThreeDimOilWaterImpes
from CellContainer import CellContainer
from Flow import Flow
from SolverSlau import SolverSlau

solver_slau = SolverSlau()
impes = ThreeDimOilWaterImpes(solver_slau)

cell_container = CellContainer()  # Проверь на счет eq_index. Внутри реализации написано чо каво
cell_container.initialize_cells()
cell = cell_container.get_cell(3, 3, 3)

cell_container.get_cell(3, 3, 3).get_cell_state_n().set_pressure_oil(200*101325)
cell_container.get_cell(3, 3, 3).get_cell_state_n_plus().set_pressure_oil(200*101325 + 1000)

cell_container.initialize_flows()
Flow.initialize_flow(cell_container)

time = impes.tau  # Сразу обозначим это как первый шаг по времени, потому что нулевой у нас есть
counter = 1
counter_write = []

while time < impes.time_max:
    delta_k = impes.generate_delta_k()
    if counter == 1:
        pass
    else:
        cell_container.equate_cell_states()

    while impes.check_norm(delta_k):
        impes.solver_slau.set_zero()
        impes.recount_properties(cell_container)
        impes.count_cells_flows(cell_container)
        impes.generate_matrix(cell_container)
        impes.solve_slau()
        delta_k = impes.solver_slau.get_result()
        impes.solver_slau.clear_result()
        impes.update_pressure(cell_container, delta_k)

    if counter % 100 == 0:
        f = open(f'govno{counter}.txt', 'w')
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_x):
                    cell = cell_container.get_cell(k, i, j)
                    f.write(str(cell.get_cell_state_n_plus().get_pressure_oil()) + '\n')
        f.close()

    impes.recount_properties(cell_container)
    impes.count_cells_flows(cell_container)
    impes.update_saturation(cell_container)
    impes.update_pressure_cap(cell_container)

    print(f'CURRENT TIME: {time / (impes.tau)} DAYS')
    time += impes.tau
    counter += 1
