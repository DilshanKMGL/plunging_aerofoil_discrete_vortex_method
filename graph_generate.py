import matplotlib.pyplot as plt
import numpy as np
import function_file as ff

file1 = open('new file.txt', 'r')
line = file1.readlines()

center = complex(line[1][10:len(line[1]) - 2])
a = float(line[2][9:])
time_step = float(line[6][21:])
pl_amp = float(line[9][11:])
pl_fre = float(line[10][11:])
iterate_number = int(line[7][21:])

# ---------------------- airfoil --------------------------
r = np.sqrt(center.imag ** 2 + (a - center.real) ** 2)
theta = np.linspace(0, 2.0 * np.pi, 20)  # 200000
cir = center + r * np.exp(1j * theta)
cir = np.array(cir)
jou = cir + a ** 2 / cir

# ---------------------- plot
vortex_iteration = 18
vortex_line = 21

for i in range(iterate_number):
    iteration = int(line[vortex_iteration][10:])
    print("iteration ", iteration)
    func = ff.plunging(pl_amp, pl_fre, time_step * (iteration - 1))
    jou_new = jou + complex(str(func[0]).replace('*', '').replace('I', 'j'))

    vor_list = line[vortex_line][1:len(line[vortex_line]) - 2].split(',')
    vortex = [complex(index.replace(' ', '').replace('I', 'j').replace('*', '')) for index in vor_list]
    real_part = [index.real for index in vortex]
    imag_part = [index.imag for index in vortex]
    plt.axis('off')
    plt.xlim(-5, 30)
    plt.ylim(-10, 10)

    plt.grid(False)
    plt.plot(jou_new.real, jou_new.imag, color='b')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.scatter(real_part, imag_part, s=2, color='g')
    heading = str(iteration)
    plt.savefig(heading)
    plt.close()

    vortex_iteration += 5
    vortex_line += 5

'''
vortex_iteration = 18
iteration = line[vortex_iteration][10:]
func = ff.plunging(pl_amp, pl_fre, time_step * (vortex_iteration - 1))
jou_new = jou + complex(str(func[0]).replace('*', '').replace('I', 'j'))

vortex_line = 27
vor_list = line[vortex_line][1:len(line[vortex_line]) - 2].split(',')
vortex = [complex(index.replace(' ', '').replace('I', 'j').replace('*', '')) for index in vor_list]
real_part = [index.real for index in vortex]
imag_part = [index.imag for index in vortex]
plt.axis('off')
plt.xlim(-5, 30)
plt.ylim(-2, 4)

plt.grid(False)
plt.plot(jou_new.real, jou_new.imag, color='b')
plt.gca().set_aspect('equal', adjustable='box')
# plt.scatter(real_part, imag_part, s=10, color='g')
plt.show()
vortex_iteration += 5
vortex_line += 5
'''
