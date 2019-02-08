from CellContainer import CellContainer
from Flow import Flow
from Layer import Layer
from SolverSlau import SolverSlau
from ThreeDimOilWaterImpes import ThreeDimOilWaterImpes
from prettytable import PrettyTable

solver_slau = SolverSlau()
impes = ThreeDimOilWaterImpes(solver_slau)

cell_container = CellContainer()  # Проверь на счет eq_index. Внутри реализации написано чо каво
cell_container.initialize_cells()

"""
cell = cell_container.get_cell(2, 2, 2)
cell.get_cell_state_n().set_pressure_water(200 * 101325)
cell.get_cell_state_n_plus().set_pressure_water(200 * 101325 + 1000)

cell.get_cell_state_n().set_pressure_oil(200 * 101325)
cell.get_cell_state_n_plus().set_pressure_oil(200 * 101325 + 1000)
"""

Flow.initialize_flow(cell_container)

time = impes.tau  # Сразу обозначим это как первый шаг по времени, потому что нулевой у нас есть
counter = 1
counter_write = []

t_debit = PrettyTable(['Time, days', 'Q_mass_OIL, kg/sec', 'Q_mass_WATER, kg/sec', 'Q_vol_OIL, cub.met/sec', 'Q_vol_WATER, cub.met/sec'])

while time < impes.time_max:
    delta_k = impes.generate_delta_k()
    if counter == 1:
        pass
    else:
        cell_container.equate_cell_states()
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_y):
                    cell = cell_container.get_cell(k, i, j)
                    cell.get_cell_state_n_plus().set_pressure_water(cell.get_cell_state_n().get_pressure_water() + 1000.0)
                    cell.get_cell_state_n_plus().set_pressure_oil(cell.get_cell_state_n().get_pressure_oil() + 1000.0)

    while impes.check_norm(delta_k):
        impes.solver_slau.set_zero()
        impes.recount_properties(cell_container)
        impes.count_cells_flows(cell_container)
        impes.generate_matrix(cell_container)
        impes.solve_slau()
        delta_k = impes.solver_slau.get_result()
        impes.solver_slau.clear_result()
        """
        if impes.check_pressure_convergence(cell_container, delta_k):
            impes.tau = impes.tau / 2.0

            for k in range(Layer.N_z):
                for i in range(Layer.N_x):
                    for j in range(Layer.N_y):
                        cell = cell_container.get_cell(k, i, j)
                        cell.get_cell_state_n_plus().set_equals_to(cell.get_cell_state_n())
                        cell.get_cell_state_n_plus().set_pressure_water(cell.get_cell_state_n().get_pressure_water() + 1000.0)
                        cell.get_cell_state_n_plus().set_pressure_oil(cell.get_cell_state_n().get_pressure_oil() + 1000.0)

            print("be")
        else:
        """
        print("govno")
        impes.update_pressure(cell_container, delta_k)

    if counter % 10 == 0:
        f = open(f'data_{counter}.txt', 'w')
        f_deb = open(f'debit.txt', 'w')
        t = PrettyTable(['Pressure_OIL, atm', 'Pressure_WATER, atm', 'Pressure_CAP, atm', 'Saturation_OIL', 'Saturation_WATER'])
        for k in range(Layer.N_z):
            for i in range(Layer.N_x):
                for j in range(Layer.N_x):
                    cell = cell_container.get_cell(k, i, j)
                    cell_state = cell.get_cell_state_n_plus()
                    p_oil = cell_state.get_pressure_oil()
                    p_water = cell_state.get_pressure_water()
                    p_cap = cell_state.get_pressure_cap()
                    s_oil = cell_state.get_s_oil()
                    s_water = cell_state.get_s_water()
                    t.add_row([p_oil / Layer.atm, p_water / Layer.atm, p_cap / Layer.atm, s_oil, s_water])

                    if cell.has_well:
                        well = cell.well
                        q_oil = well.get_oil_well_index() * (p_oil - well.p_well)
                        q_water = well.get_water_well_index() * (p_water - well.p_well)
                        q_vol_oil = q_oil / cell_state.get_ro_oil()
                        q_vol_water = q_water / cell_state.get_ro_water()
                        cur_time = time / impes.tau_default
                        t_debit.add_row([cur_time, q_oil, q_water, q_vol_oil, q_vol_water])

        f_deb.write(str(t_debit))
        f.write(str(t))
        f.close()
        f_deb.close()

    impes.recount_properties(cell_container)
    impes.count_cells_flows(cell_container)
    impes.update_saturation(cell_container)


    impes.update_pressure_cap(cell_container)

    print(f'CURRENT TIME: {time / impes.tau_default} DAYS')
    time += impes.tau
    impes.tau = min(impes.tau * 2.0, impes.tau_default)
    print(f'TIMESTEP: {impes.tau}')
    counter += 1

