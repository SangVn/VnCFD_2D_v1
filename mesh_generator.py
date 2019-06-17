# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

import  numpy as np

# Tạo lưới bài toán Mach 2
def generate_mesh_M2(Nj, Ni):
    # Kích thước vùng tính toán
    ly, lx = 4.0, 10.0
    
    # Tạo mảng 3 chiều để chứa tọa độ các điểm lưới 
    points = np.zeros((Nj, Ni, 2))
    
    # tọa độ x tại các điểm lưới
    dx = lx / Ni
    x = np.linspace(0, lx, Ni)
    
    # tọa độ y của biên dưới
    y0 = np.zeros(Ni)

    # index i tương ứng vị trí x = 2, 4 trên biên dưới
    i2 = int(2./dx)
    i4 = int(4./dx)

    y0[i2:i4] = (x[i2:i4]-2.)*np.tan(np.pi/12)
    y0[i4:] = 2.0*np.tan(np.pi/12)

    # khoảng cách dy giữa hai điểm cùng cột
    dy = np.array([(ly-y)/(Nj-1) for y in y0])
    
    # xác định tọa độ (x, y) của từng điểm 
    for j in range(Nj):
        for i in range(Ni):
            points[j, i, 0] = x[i]
            points[j, i, 1] = y0[i]+j*dy[i]

    return points


# Tạo lưới bài toán dòng chảy quanh hình trụ
def generate_mesh_cylinder(Nj, Ni):
    # chia góc 2Pi ra thành Ni điểm lưới
    alpha = np.linspace(0.0, -2 * np.pi, Ni)

    # bán kính tại các điểm lưới
    r = np.zeros(Nj)
    r[0] = 1.0  # bán kính hình trụ
    dr = 1e-2  # kích thước ('dộ dày') ô lưới đầu tiên sát biên (bề mặt trụ)
    ratio = 1.2  # tỷ lệ tăng kích thước ô lưới
    for j in range(1, Nj):
        r[j] = r[j - 1] + dr
        dr *= ratio

        # tạo độ điểm lưới
    points = np.zeros((Nj, Ni, 2))
    for j in range(Nj):
        for i in range(Ni):
            points[j, i, 0] = r[j] * np.cos(alpha[i])
            points[j, i, 1] = r[j] * np.sin(alpha[i])
    return points