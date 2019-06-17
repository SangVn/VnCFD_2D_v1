# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from lib.data import Cells, Sides
from setup import init_field, mesh_filename, field_filename
from lib.solver import eu_solver
from lib.functions import import_mesh, read_field, write_cell_data, write_block_data, show_field, show_point_data

def run():
    # nhập lưới
    Nj, Ni, points = import_mesh(mesh_filename)
    # Nj, Ni, points = import_mesh(mesh_filename, dl=',')

    # tạo dữ liệu cells, sides
    cells = Cells(Nj-1, Ni-1, points)
    sides = Sides(points, cells)

    # thiết lập trường ban đầu
    iter, time = 0, 0.0
    init_field(cells)
    write_cell_data(cells, iter, time, field_filename)

    # đọc trường ban đầu
    iter, time = read_field(cells, field_filename)

    # gọi eu_solver
    eu_solver(cells, sides, iter, time)
    write_block_data(cells, points, iter, time, field_filename.replace('cell_data.dat', 'block_data.dat'))

    # biểu diễn kết quả, lưu ảnh
    show_field(cells, points)


if __name__ == '__main__':
    run()
    # show_point_data(41, 101, 'data/mach_2_point_data.txt')
    # show_point_data(41, 61, 'data/cylinder_point_data.txt')
    # show_point_data(47, 103, 'data/naca0012_point_data.txt')
    # show_point_data(51, 157, 'data/CrewDragon_point_data.txt')