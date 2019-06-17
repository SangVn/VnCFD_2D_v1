# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from numpy import zeros, array, loadtxt
import matplotlib.pyplot as plt

# Hàm xuất lưới
def export_mesh(Nj, Ni, points, file_name):
    f = open(file_name, 'w')
    f.write('TITLE = "vncfd python"\n')
    f.write('VARIABLES = "X", "Y"\n')
    f.write('ZONE T="1", I= %d, J= %d\n' % (Ni, Nj))
    for j in range(Nj):
        for i in range(Ni):
            f.write('%f %f\n' % (points[j, i, 0], points[j, i, 1]))
    f.close()


# Hàm đọc lưới
def import_mesh(file_name, dl=' '):
    print('\nImport mesh from: %s\n' % file_name)
    f = open(file_name, 'r')

    # đọc và hiện ra màn hình 3 dòng đầu
    for i in range(3):
        line = f.readline()
        print(line)

    # lấy giá trị Ni, Nj
    words = line.split()  # chia dòng cuối ra thành các từ riêng biệt bằng dấu cách ' '
    Nj = int(words[-1])  # từ cuối cùng là Nj
    Ni = int(words[-3].replace(',', ''))  # từ thứ 3 tứ cuối lên bỏ dấu ',' là Ni

    f.close()

    # đọc tọa độ các điểm lưới bằng hàm loadtxt, bỏ 3 hàng đầu
    # dùng reshape để chuyển mảng về 3 chiều
    points = loadtxt(file_name, skiprows=3, usecols=(0,1), delimiter=dl).reshape((Nj, Ni, 2))

    return Nj, Ni, points

# Hàm P2U: xác định U từ P
gamma = 1.4 # số mũ đoạn nhiệt

def P2U(P):
    U = zeros(4)
    U[0] = P[0]
    U[1] = P[0] * P[1]
    U[2] = P[0] * P[2]
    U[3] = P[3] / (gamma - 1) + 0.5 * P[0] * (P[1] ** 2 + P[2] ** 2)
    return U

# Hàm P2F: tính dòng qua mặt (công thức ở bài 18)
def P2F(P, side):
    n = side.normal    #vector pháp tuyến đơn vị của mặt
    vn = n.dot(P[1:3]) #vận tốc vuông góc bề mặt V.n
    F = zeros(4)
    F[0] = P[0] * vn
    F[1] = F[0] * P[1] + P[3] * n[0]
    F[2] = F[0] * P[2] + P[3] * n[1]
    F[3] = F[0] * (P[3] / P[0] * gamma / (gamma - 1.0) + 0.5 * (P[1] ** 2 + P[2] ** 2))
    return F * side.area

# Hàm U2P: xác định biến biên thủy P từ biến bảo toàn U
def U2P(U, P):
    P[0] = U[0]
    P[1] = U[1] / U[0]
    P[2] = U[2] / U[0]
    P[3] = (U[3] - 0.5 * P[0] * (P[1] ** 2 + P[2] ** 2)) * (gamma - 1)

# hàm tính số mach
def Mach(P):
    a = (gamma*P[3]/P[0]) ** 0.5
    u = (P[1]*P[1] + P[2]*P[2]) ** 0.5
    M = u/a
    return M

R_gas = 287.052873836 # Hàng số chất khí
# hàm xác định khối lượng riêng phụ thuộc nhiệt độ và áp suất theo phương trình trạng thái
def rho(T, p):
    rho = p/(R_gas*T)
    return rho

def Temperature(rho, p):
    T = p/(R_gas*rho)
    return T

# Lưu tọa độ tâm ô lưới và thông số dòng chảy tại đó
# Two-Dimensional Field Plots
def write_cell_data(cells, iter, time, filename):
    print('\nWrite cell data to: %s\n' % filename)
    f = open(filename, 'w')
    f.write('TITLE = "vncfd field: iter= %d, time= %f"\n' % (iter, time))
    f.write('VARIABLES = "X", "Y", "rho", "u", "v", "p", "Mach", "T"\n')
    f.write('ZONE T="1", I=%d, J=%d, DATAPACKING=POINT\n' % (cells.size[1], cells.size[0]))
    for cell in cells:
        C = cell.center
        P = cell.P
        M = Mach(P) # hàm xác định số mach
        T = Temperature(P[0], P[3])
        f.write('%f %f %f %f %f %f %f %f\n' % (C[0], C[1], P[0], P[1], P[2], P[3], M, T))
    f.close()

# Lưu tạo độ điểm lưới và thông số dòng chảy tại tâm ô lưới
# Cell-Centered Data
def write_block_data(cells, points, iter, time, filename):
    print('\nWrite block data to: %s\n' % filename)
    f = open(filename, 'w')
    f.write('TITLE = "vncfd field: iter= %d, time= %f"\n' % (iter, time))
    f.write('VARIABLES = "X", "Y", "rho", "u", "v", "p", "Mach", "T"\n')
    f.write('ZONE T="1", I=%d, J=%d, DATAPACKING=BLOCK, VARLOCATION=([3,4,5,6,7,8]=CELLCENTERED)\n' % (points.shape[1], points.shape[0]))

    X_p, Y_p = points[:, :, 0].ravel(), points[:, :, 1].ravel()
    for x in X_p: f.write('%f ' % x)
    f.write('\n')
    for y in Y_p: f.write('%f ' % y)
    f.write('\n')

    for i in range(4):
        for cell in cells: f.write('%f ' % cell.P[i])
        f.write('\n')
    for cell in cells:
        M = Mach(cell.P)
        f.write('%f ' % M)
    for cell in cells:
        T = Temperature(cell.P[0], cell.P[3])
        f.write('%f ' % T)
    f.write('\n')
    f.close()

# đọc các thông số ban đầu từ file
def read_field(cells, filename):
    print('\nRead cell data from: %s\n' % filename)

    f = open(filename, 'r')
    line = f.readline()   # đọc dòng đầu tiên chứa iter và time
    words = line.split()  # chia dòng cuối ra thành các từ riêng biệt bằng dấu cách ' '
    time = float(words[-1].replace('"', '')) # từ cuối cùng bỏ dấu '"'là time
    iter = int(words[-3].replace(',', ''))  # từ thứ 3 tứ cuối lên bỏ dấu ',' là iter
    f.close()

    data = loadtxt(filename, skiprows=3, usecols=(2, 3, 4, 5)) # rho u v p
    for i in range(cells.len):
        cells[i].P = data[i]
        cells[i].U = P2U(data[i]) # hàm tính U từ P
    return iter, time

# Hàm biểu diễn kết quả bằng pyplot: pcolor hay contourf
def show_field(cells, points, porc=0): #porc=0 - pcolor, porc=1 - contourf
    Nj, Ni = cells.size[0], cells.size[1]

    if porc == 0:
        # Nếu vẽ pcolor, pcolormesh... cần tọa độ tại đỉnh, giá trị tại tâm
        X_p, Y_p = points[:, :, 0], points[:, :, 1]
    else:
        # Nếu vẽ contour, contourf, quiver, streamplot... cần tọa độ tại tâm, giá trị tại tâm
        centers = array([cell.center for cell in cells]).reshape((Nj, Ni, 2))
        X_c, Y_c = centers[:, :, 0], centers[:, :, 1]

    fig, axs = plt.subplots(2, 2)
    titles = ['rho', 'u', 'v', 'p']
    i = 0
    for ax in axs.flat:
        value_c = array([cell.P[i] for cell in cells]).reshape((Nj, Ni))
        if porc == 0: c = ax.pcolor(X_p, Y_p, value_c)
        else: c = ax.contourf(X_c, Y_c, value_c)
        ax.set_title(titles[i])
        fig.colorbar(c, ax=ax)
        i += 1

    fig.tight_layout()
    # plt.savefig('test.png')
    plt.show()


# Chuyển dữ liệu từ dạng block_data sang point_data:
# Bước 1. dùng ParaView mở file block_data
# Bước 2. Apply filter: CellDataToPointData
# Bước 3. xuất dữ liệu tại điểm lưới: save data -> 'filename.txt'
def show_point_data(Nj, Ni, filename):
    data = loadtxt(filename, skiprows=1, delimiter=',')  # file thu được từ paraview chuyển dạng 2 sang 3
    # dòng đầu tiên: "rho", "u", "v", "p", "Mach", "T", "Points:0", "Points:1", "Points:2"
    # tương ứng các cột: 0, 1, 2, 3, 4, 5, 6, 7, 8

    x = data[:, 6].reshape((Nj, Ni)) # cột thứ 5
    y = data[:, 7].reshape((Nj, Ni)) # cột thứ 6

    fig, axs = plt.subplots(2, 2)
    titles = ['rho', 'u', 'v', 'p', 'Mach', 'T']
    iv = [0, 3, 4, 5] # chỉ số các cột muốn biểu diễn

    i = 0
    for ax in axs.flat:
        value = data[:, iv[i]].reshape((Nj, Ni))
        c = ax.contourf(x, y, value)
        ax.set_title(titles[iv[i]])
        # ax.set_xlim([-25, 15])
        # ax.set_ylim([-20, 20])
        # ax.set_xlim([-0.5, 1.5])
        # ax.set_ylim([-0.5, 0.5])
        fig.colorbar(c, ax=ax)
        i += 1

    fig.tight_layout()
    img_filename = filename.replace('.txt', '.png')
    plt.savefig(img_filename)
    plt.show()
