import numpy as np
import math


def Coil_resistance(material, l, d, nu):
    """
    Calculates the coil resistance
    ---------------
    @param material: Contour material
    @param l: Contour length
    @param d: Diameter of the conductor cross section
    @param nu: The frequency of the current in the contour
    @return: Coil resistance
    """

    Epsilon0 = 8.85419e-12
    omega = 2 * np.pi * nu
    c = 299792458

    ro = {
        'Silver': 0.016,
        'Copper': 0.0168,
        'Gold': 0.024,
        'Aluminum': 0.028,
        'Tungsten': 0.055}

    ro[material] = ro[material] * 1e-6

    delta = c * np.sqrt((2 * Epsilon0 * ro[material]) / omega)

    S_eff = np.pi * (d / 2) ** 2 - np.pi * (d / 2 - delta) ** 2

    R = (ro[material] * (l / S_eff))

    return R


def resistance_contour(l, material, d, nu):
    """
    Calculates the contour resistance
    ---------------
    @param l: Array of lengths of parallel connected coils
    @param material: Contour material
    @param d: Diameter of the conductor cross section
    @param nu: The frequency of the current in the contour
    @return: Contour resistance
    """
    Resistance_coil = []
    for i in l:
        R = Coil_resistance(material, i, d, nu)
        Resistance_coil.append(R)
    return (np.sum(list(map(lambda x: x ** (-1), Resistance_coil)))) ** (-1)


def length_circular_coils(coils):
    """
    Calculates the length of each sequentially connected turn in a circular contour
    ---------------
    @param coils: Multidimensional array with the radii of the turns
    @return l: Array with lengths of sequentially connected turns
    """
    l = []
    for coil in coils:
        l.append(2 * np.pi * np.sum(coil))
    return l


def length_square_coils(coils):
    """
    Calculates the length of each sequentially connected turn in a rectangular contour
    ---------------
    @param coils: Multidimensional array with the lengths of the sides of the turns
    @return l: Array with lengths of sequentially connected turns
    """
    l = []
    for coil in coils:
        lenght_coil = 0
        for i in coil:
            lenght_coil += 2 * np.sum(i)
        l.append(lenght_coil)
    return l


def length_piecewise_linear_coils(coils):
    """
    Calculates the length of each sequentially connected turn in a piecewise linear contour
    ---------------
    @param coils: Multidimensional array with the coordinates of the turns
    @return l: Array with lengths of sequentially connected turns
    """
    l = []
    for coil in coils:
        lenght_coil = 0
        for j, turn in enumerate(coil):
            for i in range(len(turn)):
                try:
                    lenght_coil += np.sqrt(
                        (coil[j][i][0] - coil[j][i + 1][0]) ** 2 + (coil[j][i][1] - coil[j][i + 1][1]) ** 2)
                except IndexError:
                    lenght_coil += np.sqrt((coil[j][0][0] - coil[j][i][0]) ** 2 + (coil[j][0][1] - coil[j][i][1]) ** 2)
        l.append(lenght_coil)
    return l


def split(turns, freq):
    """
    Function that breaks up the original array of turn radiuses into a 2-dimensional array of multiple separate coils,
    that are to be connected as parallel circuit. This is required to match the maximum possible length of the wires
    which is for magnetostatic approximation to work.
    @param turns: initial array of radiuses
    @param freq: frequency of the EM wave, is used to determine the maximum possible wire length
    @return: the broken up array of radiuses
    """
    # This code works best if the :turns: array is sorted, which it is.
    # freq = freq * 1_000_000
    length_of_wave = c / freq
    max_length = length_of_wave / 6
    n = math.ceil((sum(turns) * 2 * math.pi) / max_length)

    res = []
    for i in range(n):
        res.append([])

    for turn in turns:
        res[0].append(turn)
        res = sorted(res, key=lambda x: sum(x))

    return res
