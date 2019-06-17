# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from numpy import array
from lib.fluxes import flux_roe
from lib.boco import boco_dict, set_freestream
from lib.functions import P2U, rho

# cylinder, naca
P_freestream = array([rho(293.15, 101325.0), 40, 0.0, 101325.0]) # freestream

# CrewDragon
# P_freestream = array([rho(293.15, 101325.0), 500, 0.0, 101325.0]) # freestream


def init_field(cells):
    P_t0 = P_freestream
    U_t0 = P2U(P_t0)
    for cell in cells[:]:
        cell.P = P_t0.copy()
        cell.U = U_t0.copy()
    
# điều kiện biên
def set_boco(sides):
    set_freestream(P_freestream)
    boco_name = ['joint', 'joint', 'wall', 'farfield']
    boco_func = [boco_dict[name] for name in boco_name]
    sides.set_joint(sides.boundary_0, sides.boundary_1)
    sides.set_boco(boco_func)

# các thông số khác
CFL = 0.85
time_target = None
iter_target = 10

write_field_iter = 1000 # thời điểm ghi kết quả giữa chừng
flux_func = flux_roe

mesh_filename = 'data/cylinder_mesh.dat'
field_filename = 'data/cylinder_cell_data.dat'

# mesh_filename = 'data/naca0012_mesh.dat'
# field_filename = 'data/naca0012_cell_data.dat'

# mesh_filename = 'data/CrewDragon_mesh.dat'
# field_filename = 'data/CrewDragon_cell_data.dat'