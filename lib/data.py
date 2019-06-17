# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

import numpy as np
from functions import U2P

# tâm ô lưới bằng trung bình cộng tọa độ bốn đỉnh lưới
def center(vertices):
    return sum(vertices) / 4.

# diện tích ô lưới bằng một nửa độ lớn tích có hướng hai vector đường chéo
# ~ tương đương thể tích ô lưới trong trường hợp 3 chiều
def volume(vertices):
    return abs(np.cross(vertices[0] - vertices[2], vertices[1] - vertices[3]) / 2.)

def cell_size(vertices):
    dx_vec = (vertices[1]-vertices[0] + vertices[2] - vertices[3])/2.
    dy_vec = (vertices[2]-vertices[1] + vertices[3] - vertices[0])/2.
    dx = dx_vec.dot(dx_vec)**0.5
    dy = dy_vec.dot(dy_vec)**0.5
    return min(dx, dy)

# lớp dữ liệu ô lưới:
# tại thời điểm khởi tạo, giá trị P, U chưa xác định (None)
class Cell:
    # Hàm khởi tạo, khai báo
    def __init__(self, vertices): #4 đỉnh lưới theo thứ tự ngược chiều KĐH
        self.center   = center(vertices)
        self.volume   = volume(vertices)
        self.size     = cell_size(vertices)
        self.P = None           #(rho, u, v, p)
        self.U = None           #(rho, rho*u, rho*v, rho*e); e= p/(rho*km1) + (u**2 + v**2)/2
        self.res = np.zeros(4)  #tổng các dòng đi qua các bề mặt: sum(F*S)
        self.dt  = 0.0          #bước thời gian
        
# Để xác định các ô lưới ta cần kích thước lưới NjxNi và tọa độ các điểm lưới points
class Cells():
    # Khởi tạo
    def __init__(self, Nj, Ni, points):
        self.size  = [Nj, Ni] # Kích thước lưới
        self.len   = Nj*Ni    # Tổng số ô lưới
        self.cells = []       # Dãy 1D chứa các ô lưới
        for j in range(Nj):
            for i in range(Ni):
                #ô lưới thứ (j,i) gồm 4 đỉnh [(j,i), (j,i+1), (j+1,i+1), (j+1,i)] (hình 1)
                vers = (points[j, i], points[j, i + 1], points[j + 1, i + 1], points[j + 1, i])
                # khởi tạo ô lưới và thêm vào dãy các ô lưới
                self.cells.append(Cell(vers))

    # Ba phương thức cơ bản để lấy các phần tử của dãy các ô lưới:
    def __getitem__(self, item):
        # Lấy một đoạn các ô lưới: cells[start:stop]
        if isinstance(item, slice):
            return self.cells.__getslice__(item.start, item.stop)
        
        # Lấy ô lưới thứ (j,i): cells[j,i]
        elif isinstance(item, tuple):
            j, i = item
            if(j < 0): j += self.size[0]
            if(i < 0): i += self.size[1]
            return self.cells[j*self.size[1]+i]
        
        # Lấy ô lưới thứ j*i: cells[j*i]
        elif isinstance(item, int):
            return self.cells[item]

    # tính bước thời gian cục bộ trong từng ô lưới
    def time_step_cell(selfs):
        for cell in selfs.cells:
            a = (1.4 * cell.P[3] / cell.P[0]) ** 0.5  # vận tốc âm thanh
            v = (cell.P[1] ** 2 + cell.P[2] ** 2) ** 0.5  # vận tốc dòng chảy
            cell.dt = cell.size/(v + a)

    # Xác định bước thời gian cho toàn vùng tính toán (toàn cục)
    def time_step_global(self, CFL): # CFL mặc định bằng 0.25
        # trước hết cần xác định bước thời gian trong từng ô lưới
        self.time_step_cell()
        # sau đó tìm bước thời gian nhỏ nhất trong các ô lưới
        dt = 1e6
        # tìm bước thời gian nhỏ nhất trong các ô lưới
        for cell in self.cells:
            dt = min(dt, cell.dt)          
        return CFL*dt
    
    # Thực hiện bước lặp: xác định U ở bước thời gian tiếp theo
    def new_U(self, dt):
        for cell in self.cells:
            # print('cell.ress ', cell.res)
            cell.U += dt/cell.volume*cell.res # công thức (1)
            cell.res[:] = 0.0                 # sau khi xác định U, đưa giá trị res về 0.0

    # Xác định P từ U
    def new_P(self):
        for cell in self.cells:
            U2P(cell.U, cell.P) # Hàm U2P

#sử dụng tích vô hướng hai vector .dot để tính chiều dài vector
#~ tương đương diện tích bề mặt trong trường hợp 3 chiều
def area(side_vec):
    return (side_vec.dot(side_vec))**0.5

#xác định vector pháp tuyến của bề mặt: S*n, với n - vector pháp tuyến đơn vị
def normal(side_vec):
    return np.array([side_vec[1], -side_vec[0]])

#định nghĩa lớp bề mặt
class Side:
    def __init__(self, side_vec):
        self.area   = area(side_vec) # diện tích bề mặt
        self.normal = normal(side_vec)/self.area   #vector pháp tuyến đơn vị
        self.cells = None   #hai ô lưới hai bên trái phải, tạm thời chưa xác định
        # Quy ước: self.cells[0] - ô bên trái, self.cells[1] - ô bên phải

'''
Quy ước:
                 boundary_3
                  <--------    
              ^             ^
    boundary_0| inner_sides |boundary_1
                  <--------
                 boundary_2  
'''

class Sides:   
    # để khởi tạo Sides cần tọa độ điểm lưới về địa chỉ các ô lưới
    def __init__(self, points, cells):
        #xác định các vector bề mặt từ các điểm lưới 
        sides_i = points[:, :-1] - points[:, 1:]   # sides nằm ngang
        sides_j = points[1:] - points[:-1]         # sides thẳng đứng

        # Biên_0 gồm các mặt ở cột đầu sides_j
        # Ô lưới bên trái không xác định, bên phải là các ô ở cột đầu tiên 
        self.boundary_0 = []
        for j in range(sides_j.shape[0]):
            side = Side(sides_j[j, 0])
            side.cells = [None, cells[j, 0]]
            self.boundary_0.append(side)
            # print('B0', side.normal)

        # Biên_1 gồm các mặt ở cột cuối sides_j
        # Ô lưới bên phải không xác định, bên trái là các ô ở cột cuối
        self.boundary_1 = []
        for j in range(sides_j.shape[0]):
            side = Side(sides_j[j, -1])
            side.cells = [cells[j, -1], None]
            self.boundary_1.append(side)
            # print('B1', side.normal)

        # Biên_2 gồm các mặt ở hàng đầu sides_i
        # Ô lưới bên phải không xác định, bên trái là các ô ở hàng đầu
        self.boundary_2 = []
        for i in range(sides_i.shape[1]):
            side = Side(sides_i[0, i])
            side.cells = [None, cells[0, i]]
            self.boundary_2.append(side)
            # print('B2', side.normal)

        # Biên_3 gồm các mặt ở hàng cuối sides_i
        # Ô lưới bên trái không xác định, bên phải là các ô ở hàng cuối
        self.boundary_3 = []
        for i in range(sides_i.shape[1]):
            side = Side(sides_i[-1, i])
            side.cells = [cells[-1, i], None]
            self.boundary_3.append(side)
            # print('B3', side.normal)

        # Những hàng còn lại bên trong sides_i
        self.inner_sides = []
        for j in range(1, sides_i.shape[0] - 1):
            for i in range(sides_i.shape[1]):
                side = Side(sides_i[j, i])
                side.cells = [cells[j - 1, i], cells[j, i]]
                self.inner_sides.append(side)

        # Những cột còn lại bên trong sides_j
        for i in range(1, sides_j.shape[1] - 1):
            for j in range(sides_j.shape[0]):
                side = Side(sides_j[j, i])
                side.cells = [cells[j, i - 1], cells[j, i]]
                self.inner_sides.append(side)
    
    # Hàm đặt điều kiện biên 
    def set_boco(self, boco):
        self.boco = boco # boco là tập hợp các điều kiện biên

    # Trường hợp điều kiện biên tuần hoàn
    def set_joint(self, boundary_0, boundary_1):
        for side_0, side_1 in zip(boundary_0, boundary_1):
            side_0.cells[0] = side_1.cells[0]
            self.inner_sides.append(side_0)

    # Hàm xác định dòng qua biên
    def flux_boundaries(self):
        # Tương ứng 4 biên, có 4 điều kiện biên (trường hợp đơn giản)
        # Các chỉ số 1, 0, 1, 0 là chỉ số của ô lưới kề biên side.cells[id]
        self.boco[0](self.boundary_0, 1) # ô lưới kề biên ở bên phải
        self.boco[1](self.boundary_1, 0) # ô lưới kề biên ở bên trái
        self.boco[2](self.boundary_2, 1) # ô lưới kề biên ở bên phải
        self.boco[3](self.boundary_3, 0) # ô lưới kề biên ở bên trái

    # Hàm xác định dòng qua các biên bên trong
    def flux_inner_sides(self, flux_func): #flux_func: hàm tính dòng qua từng mặt (side)
        for side in self.inner_sides:
            F = flux_func(side, side.cells[0].P, side.cells[1].P)
            # Như đã quy ước ở trên:
            side.cells[0].res -= F # ô bên trái -
            side.cells[1].res += F # ô bên phải +
            
