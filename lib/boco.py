# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from functions import P2F

# thông số dòng tự do
P_freestream = None

def set_freestream(P):
    global P_freestream
    P_freestream = P

# Hàm sign_ic: xác định `thêm hay bớt` dòng ở ô lưới kề bên side
def sign_ic(ic): #ic ~ index_cell; left_cell: ic = 0; right_cell: ic = 1
    if ic == 0: return -1.0 # res -= flux
    else: return 1.0        # res += flux

# Điều kiện dòng chảy vào
# P_side = P_freestream
def supersonic_inflow(boundary, ic):
    for side in boundary:
        F = P2F(P_freestream, side)
        side.cells[ic].res += sign_ic(ic)*F

# Điều kiện dòng chảy ra
# P_side = P_in
def supersonic_outflow(boundary, ic):
    for side in boundary:
        F = P2F(side.cells[ic].P, side)
        side.cells[ic].res += sign_ic(ic)*F

# Điều kiện biên wal_no_slip
# P_side: u=v=0, dp/dn=0, drho/dn=0
def wall_no_slip(boundary, ic):
    for side in boundary:
        P_side = [side.cells[ic].P[0], 0.0, 0.0, side.cells[ic].P[3]]
        F = P2F(P_side, side)
        side.cells[ic].res += sign_ic(ic)*F

# Điều kiện biên farfield
def farfield(boundary, ic):
    for side in boundary:
        V_in = side.cells[ic].P[1:3] # vận tốc tại ô lưới phía trong biên
        Vn = -sign_ic(ic)*side.normal.dot(V_in) # hướng vận tốc so với phương pháp tuyến (quay ra) của biên
        if Vn >= 0: # chảy ra
            F = P2F(side.cells[ic].P, side)
        else: # chảy vào
            F = P2F(P_freestream, side)
        side.cells[ic].res += sign_ic(ic) * F

# Điều kiện biên joint
def joint(boundary, ic):
    pass

# ta có `từ điển` điều kiện biên: <tên gọi> : <hàm đkb>
boco_dict = {'supersonic_inflow': supersonic_inflow, 'supersonic_outflow': supersonic_outflow, 'wall': wall_no_slip, \
             'farfield': farfield, 'joint': joint}
