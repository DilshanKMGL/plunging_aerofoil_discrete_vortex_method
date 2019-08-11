import numpy as np
import sympy as sp
import cmath


# -------------------------- general mapping function ------------------------------------------------------------------
def mapzeta(z_coordinate, a):
    """

    complex number / list :param z_coordinate: can be a complex number or complex number list
    float
    :param z_coordinate: coordinate value in z plane
    :param a: joukowski parameter
    complex number or complex number list
    :return: complex number or complex number list
    """
    z, zeta = sp.symbols('z, zeta')
    function = z - (zeta + pow(a, 2) / zeta)
    if isinstance(z_coordinate, list):
        z_coordinate = [function.subs(z, index) for index in z_coordinate]
        answer = [sp.solve(index, zeta) for index in z_coordinate]
    else:
        z_coordinate = function.subs(z, z_coordinate)
        answer = sp.solve(z_coordinate, zeta)
    return answer


def check_in_zeta(z_coordinate, center_coor, a):
    """
    :param z_coordinate:
    :param center_coor:
    :param a:
    :return:
    """
    radius = np.sqrt(center_coor.imag ** 2 + (a - center_coor.real) ** 2)
    distance = [cmath.polar(index - center_coor)[0] for index in z_coordinate]
    if distance[0] < radius: # changed that position
        return z_coordinate[0]
    else:
        return z_coordinate[1]


def plunging(amplitude, f, current_time):
    """
    :param current_time: time of the current iterating
    :param amplitude: amplitude of the plunging
    :param f: frequency of the plunging
    :return: plunging function, plunging velocity
    """
    t = sp.symbols('t')
    fun = 1j * amplitude * sp.sin(2 * sp.pi * f * t)
    vel = sp.diff(fun, t)
    fun = fun.subs(t, current_time).evalf()
    vel = vel.subs(t, current_time).evalf()
    return [fun, vel]


def get_trailing_edge(center_coor, a, plung_pos):
    """
    This should be revised as the plunging motion is going to be modeled
    :param center_coor: center coordinate of joukowski airfoil
    :param plung_pos: current position in plunging
    :param a: joukowski parameter
    :return: z coordinate, zeta coordinate
    """
    z = 2 * a + 1j * sp.im(plung_pos)
    zeta = mapzeta(z, a)
    if plung_pos == 0:
        return [z, zeta[0]]
    else:
        zeta = check_in_zeta(zeta, center_coor, a)
        return [z, zeta]


# --------------------------- vortex calculation -----------------------------------------------------------------------
def new_vortex_position(trailing_edge, center, angle, distance, a):
    """
    :param center: center of airfoil
    :param a: joukowski parameter
    :param trailing_edge: current trailing edge of the airfoil
    :param angle: vortex position angle
    :param distance: vortex position distance relative to the chord length
    :return: z coordinate zeta coordinate
    """
    z = trailing_edge + distance * sp.exp(-1j * sp.rad(angle)).evalf()
    zeta = mapzeta(z, a)
    zeta = check_in_zeta(zeta, center, a)
    return [z, zeta]


def create_circulation(circle_center):
    """
    :param circle_center: center of the airfoil
    :return: circulation function
    get the circulation function around airfoil
    """
    zeta, vor = sp.symbols('zeta, vor', real=False)
    func = -1j * vor * sp.log(zeta - circle_center) / (2 * sp.pi)
    return func


def get_circulation(circle_center, strength):
    zeta = sp.symbols('zeta', real=False)
    func = -1j * strength * sp.log(zeta - circle_center) / (2 * sp.pi)
    return func


def create_vortex(circle_center, vortex_center, sum_strength, r):
    """

    :param circle_center: center of the airfoil
    :param vortex_center: center of the vortex
    :param sum_strength: summation of the currently emit vortex strength
    :param r: radius of circle
    :return: function for new shed vortex and image
    """
    zeta, vor = sp.symbols('zeta, vor', real=False)
    func1 = -1j * (- sum_strength - vor) * sp.log(zeta - vortex_center) / (2 * sp.pi)
    func2 = 1j * (- sum_strength - vor) * sp.log((r ** 2 / (zeta - circle_center)) - vortex_center) / (2 * sp.pi)
    func = func1 + func2
    return func


def get_vortex(circle_center, vortex_center, r, strength):
    """
    :param circle_center: center of the airfoil
    :param vortex_center: center of the vortex
    :param r: radius of the vortex
    :param strength: strength of the vortex
    :return: function of vortex with image
    """
    zeta = sp.symbols('zeta', real=False)
    func1 = -1j * strength * sp.log(zeta - vortex_center) / (2 * sp.pi)
    func2 = 1j * strength * sp.log((r ** 2 / (zeta - circle_center)) - vortex_center) / (2 * sp.pi)
    func = func1 + func2
    return func


# --------------------------- freestream calculation -------------------------------------------------------------------
def get_freestream(circle_center, r, vel, aoa):
    """
    :param circle_center: center of the airfoil
    :param r: radius of the airfoil
    :param vel: freestream velocity
    :param aoa: angle of attack
    :return:
    """
    zeta = sp.symbols('zeta', real=False)
    func1 = vel * zeta * (sp.cos(aoa) - 1j * sp.sin(aoa))
    func2 = vel * (r ** 2 / (zeta - circle_center)) * (sp.cos(aoa) + 1j * sp.sin(aoa))
    return func1 + func2


# --------------------------- equation solving -------------------------------------------------------------------------
def calculate_circulation(func, tev_edge):
    """

    :param func: function to solve
    :param tev_edge: trailing edge vortex coordinates
    :return: circulation around the airfoil
    """
    zeta, vor = sp.symbols('zeta, vor', real=False)
    func = sp.diff(func, zeta)
    func = func.subs([(zeta, tev_edge)])
    answer = sp.solve(func)
    return answer[0]


# ---------------------------- writing a file
def make_file(a, center, velocity, aoa, t_step, iteration, distance, angle, pl_amp, pl_f):
    heading = 'new file.txt'
    file1 = open(heading, 'w')

    file1.write('joukowski parameters\n')
    file1.write('center - ' + str(center) + '\n')
    file1.write('a      - ' + str(a) + '\n')

    file1.write('Flow field\n')
    file1.write('velocity           - ' + str(velocity) + '\n')
    file1.write('angle of attack    - ' + str(aoa) + '\n')
    file1.write('time step          - ' + str(t_step) + '\n')
    file1.write('iteration          - ' + str(iteration) + '\n')

    file1.write('plunging parameters\n')
    file1.write('amplitude - ' + str(pl_amp) + '\n')
    file1.write('frequency - ' + str(pl_f) + '\n')

    file1.write('new vortex\n')
    file1.write('distance   - ' + str(distance) + '\n')
    file1.write('angle      - ' + str(angle) + '\n')

    file1.write('circulation list\n')
    file1.write('vortex strength list\n')
    file1.write('vortex z coordinate list\n')
    file1.write('vortex zeta coordinate list\n')


def write_array(circulation, vortex_strength, vortex_z, vortex_zeta, i):
    heading = 'new file.txt'
    file1 = open(heading, "a+")
    file1.write('Iteration ' + str(i + 1) + '\n')
    file1.write(str(circulation))
    file1.write('\n')
    file1.write(str(vortex_strength))
    file1.write('\n')
    file1.write(str(vortex_z))
    file1.write('\n')
    file1.write(str(vortex_zeta))
    file1.write('\n')
    file1.close()
