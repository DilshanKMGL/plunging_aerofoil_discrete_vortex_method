import sympy as sp
import numpy as np
import time
import function_file as ff
import cmath

# -------------- parameters ------------------
#
# --------------------------------------------
start = time.time()
# plunging parameters
pl_amplitude = 1  # 5
pl_frequency = 5  # 0.05

# joukowski paramters
center = -0.2 + 1j * 0.2
joukowski_parameter = 2
radius = np.sqrt(center.imag ** 2 + (joukowski_parameter - center.real) ** 2)
free_velocity = 20.0
aoa = 2.0  # degree
aoa = sp.rad(aoa).evalf()

# time
time_step = 0.005
current_time = 0.00
iteration = 100

# new vortex
distance = 0.015  # portion of the chord length
angle = 0  # degrees

# data store
circulation = []
vortex_strength = []
vortex_z = []
vortex_zeta = []

# derivatives
zeta = sp.symbols('zeta', real=False)
func = zeta + pow(joukowski_parameter, 2) / zeta
first_derive = 1 / sp.diff(func, zeta)
second_derive = sp.diff(1 / sp.diff(func, zeta), zeta) / sp.diff(func, zeta)

# --------------------------------------------------------------------------------
ff.make_file(joukowski_parameter, center, free_velocity, aoa, time_step, iteration, distance, angle,
             pl_amplitude, pl_frequency)
# --------------------------------------------------------------------------------
for iterate in range(iteration):
    print('Iteration - ' + str(iterate + 1))
    # --------------------------------------------------------------------------------
    # aoa should calculate before any function run
    # --------------------------------------------------------------------------------

    fun_1 = ff.plunging(pl_amplitude, pl_frequency, current_time)  # return position and velocity
    func = cmath.polar(fun_1[1] + free_velocity)  # aoa should calculate before any function run
    velocity = func[0]
    aoa = func[1]
    
    z_coor, zeta_coor = ff.get_trailing_edge(center, joukowski_parameter, fun_1[0])  # return trailing edge - z, zeta)
    answer = ff.new_vortex_position(z_coor, center, angle, distance, joukowski_parameter)
    vortex_z.append(answer[0])
    vortex_zeta.append(answer[1])
    # create function to make complex potential function
    vortex_function = ff.create_circulation(center)
    vortex_sum = sum(vortex_strength)
    vortex_function += ff.create_vortex(center, vortex_zeta[-1], vortex_sum, radius)
    vortex_function += ff.get_freestream(center, radius, velocity, aoa)
    temp_func = [ff.get_vortex(center, vortex_zeta[vortex], radius, vortex_strength[vortex])
                 for vortex in range(len(vortex_strength))]  # all the vortex functions summation
    vortex_function += sum(temp_func)
    new_circulation = ff.calculate_circulation(vortex_function, zeta_coor)
    circulation.append(sp.re(new_circulation))  # append only real part of the circulation
    vortex_strength.append(sp.re(-vortex_sum - new_circulation))  # append only real part of the vortex strength

    # --------------------------------------------------------------------------------
    # write in csv file
    ff.write_array(circulation, vortex_strength, vortex_z, vortex_zeta, iterate)
    # --------------------------------------------------------------------------------

    current_time += time_step

    # move vortex
    vortex_function = sum(temp_func) + ff.get_circulation(center, circulation[-1]) + \
                      ff.get_vortex(center, vortex_zeta[-1], radius, vortex_strength[-1]) + \
                      ff.get_freestream(center, radius, velocity, aoa)
    vortex_function = sp.diff(vortex_function.evalf(), zeta)
    velocity_function = [(vortex_function -
                          sp.diff(ff.get_vortex(center, vortex_zeta[index],
                                                radius, vortex_strength[index]), zeta)).evalf()
                         for index in range(len(vortex_strength))]
    velocity_function = np.multiply(velocity_function, first_derive) - (np.multiply(vortex_strength, 1j) / sp.pi *
                                                                        second_derive / first_derive)
    velocity_function = [velocity_function[index].subs(zeta, vortex_zeta[index]).evalf()
                         for index in range(len(velocity_function))]
    velocity_function = np.conj(velocity_function)
    vortex_z = vortex_z + velocity_function * time_step
    vortex_z = vortex_z.tolist()
    vortex_zeta = ff.mapzeta(vortex_z, joukowski_parameter)
    vortex_zeta = [ff.check_in_zeta(index, center, joukowski_parameter) for index in vortex_zeta]

print(time.time() - start)
