# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from functions import write_cell_data
import setup
# import setup_mach_2 as setup

# Hàm eu_solver thực hiện các bước lặp để tìm nghiệm
# Biến đầu vào gồm có: các ô lưới, các mặt, số vòng lặp, thời gian lúc ban đầu
def eu_solver(cells, sides, iter, time):
    # thiết lập điều kiện biên
    setup.set_boco(sides)

    # tính theo thời gian
    if setup.time_target is not None:
        while(time < setup.time_target):
            iter += 1  # tăng số vòng lặp lên 1 đơn vị
            # tính bước thời gian
            dt = cells.time_step_global(setup.CFL)  # có thể thiết lập dt trong setup: dt = setup.dt
            if (time + dt > setup.time_target): dt = setup.time_target - time
            time += dt
            print('iteration: %d, dt: %f, time: %f' % (iter, dt, time))
            iter, time = iteration(cells, sides, iter, time, dt)

    # tính theo số vòng lặp
    elif setup.iter_target is not None:
        while(iter < setup.iter_target):
            iter += 1  # tăng số vòng lặp lên 1 đơn vị
            # tính bước thời gian
            dt = cells.time_step_global(setup.CFL)  # có thể thiết lập dt trong setup: dt = setup.dt
            time += dt
            print('iteration: %d, dt: %f, time: %f' % (iter, dt, time))
            iter, time = iteration(cells, sides, iter, time, dt)

    # ghi lại kết quả cuối cùng
    write_cell_data(cells, iter, time, setup.field_filename)


def iteration(cells, sides, iter, time, dt):
    # tính dòng qua các bề mặt
    sides.flux_boundaries()  # dòng qua các mặt biên
    sides.flux_inner_sides(setup.flux_func)  # dòng qua các mặt trong

    # tính giá trị U mới
    cells.new_U(dt)

    # tính giá trị P mới
    cells.new_P()

    # nếu số vòng lặp bằng write_field_iter thì ghi lại kết quả
    if setup.write_field_iter is not None and not (iter % setup.write_field_iter):
        write_cell_data(cells, iter, time, setup.field_filename)

    return iter, time
