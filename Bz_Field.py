import numpy as np
import math
from scipy.special import ellipk, ellipkm1, ellipe
from utilities import Radii_in_sides_square, Radii_in_coords


def Bz_segment(start_point, end_point, I, cp, calc_radius, height):
    """
    Calculates z-component of B field of single-segment
    ---------------
    @param start_point: The beginning of the segment
    @param end_point: The end of the segment
    @param I: Current in the contour
    @param cp: Calculation domain points
    @param calc_radius: Calculation domain length
    @return: Z-component of B field of single-segment
    """
    x1, y1 = start_point[0], start_point[1]
    x2, y2 = end_point[0], end_point[1]
    mu0 = np.pi * 4e-7
    C = mu0 * I / (4 * np.pi)

    x = np.linspace(-calc_radius, calc_radius, cp)
    xv, yv = np.meshgrid(x, x)

    if x1 != x2 and y1 != y2:

        k = (y2 - y1) / (x2 - x1)
        b = y2 - k*x2

        alpha = np.sqrt(k**2 + 1)
        betta = (xv + k*yv - b*k) / (alpha**2)
        gamma = yv - (k*xv + b)
        delta = height**2 + (gamma / alpha)**2

        Bz_segment_1 = C * ((x1 - betta) * gamma) / (delta * np.sqrt((alpha * (x1 - betta))**2 + delta))
        Bz_segment_2 = C * ((x2 - betta) * gamma) / (delta * np.sqrt((alpha * (x2 - betta))**2 + delta))

        return Bz_segment_2 - Bz_segment_1

    elif x1 == x2 and y1 != y2:

        alpha = height**2 + (x1 - xv)**2

        Bz_segment_1 = C * ((x1 - xv) * (y1 - yv)) / (alpha * np.sqrt((y1 - yv)**2 + alpha))
        Bz_segment_2 = C * ((x1 - xv) * (y2 - yv)) / (alpha * np.sqrt((y2 - yv)**2 + alpha))

        return Bz_segment_2 - Bz_segment_1

    elif x1 != x2 and y1 == y2:

        alpha = height ** 2 + (y1 - yv) ** 2

        Bz_segment_1 = C * ((yv - y1) * (x1 - xv)) / (alpha * np.sqrt((x1 - xv) ** 2 + alpha))
        Bz_segment_2 = C * ((yv - y1) * (x2 - xv)) / (alpha * np.sqrt((x2 - xv) ** 2 + alpha))

        return Bz_segment_2 - Bz_segment_1

def Bz_piecewise_linear_contour_single(coords,  I, cp, calc_radius, height):
    """
    Calculates Bz field of piecewise linear coil
    ---------------
    @param coords: Coordinates of the contour corners
    @param I: Current in the contour
    @param cp: Calculation domain points
    @param calc_radius: Calculation domain length
    @return: Z-component B of the field of single coil
    """

    Bz_piecewise_linear_contour_single = np.zeros((cp, cp))
    for i in range(len(coords)):
        try:
            Bz_piecewise_linear_contour_single += Bz_segment(start_point=coords[i],
                                                             end_point=coords[i + 1],
                                                             I=I,
                                                             cp=cp,
                                                             calc_radius=calc_radius,
                                                             height=height)
        except IndexError:
            Bz_piecewise_linear_contour_single += Bz_segment(start_point=coords[i],
                                                             end_point=coords[0],
                                                             I=I,
                                                             cp=cp,
                                                             calc_radius=calc_radius,
                                                             height=height)

    return Bz_piecewise_linear_contour_single


def Bz_piecewise_linear_contour(R, coords,  I, spacing, cp, height, direction=False):
    """
    Calculates the Bz field for a piecewise linear contour
    ---------------
    @param R: Array of radii
    @param coords: Coordinates of the largest turn
    @param I: Current in the contour
    @param spacing: Spacing between coil and the calculation domain boundary
    @param cp: Calculation domain points
    @param direction: The direction of the current along the contour. If the current flows clockwise, then by default this value is True
    @return: Z-component B of the field of a piecewise linear contour
    """
    if not direction:
        I = -I

    I = np.sqrt(2)*I

    list_of_coords = Radii_in_coords(R, coords)

    l = []
    for i in range(len(coords)):
        l.append(np.sqrt((coords[i][0])**2 + (coords[i][1])**2))

    calc_radius = np.amax(l) * spacing

    Bz_piecewise_linear_contour = np.zeros((cp, cp))

    for coil in list_of_coords:
        Bz_piecewise_linear_contour += Bz_piecewise_linear_contour_single(coords=coil,
                                                                          I=I,
                                                                          cp=cp,
                                                                          calc_radius=calc_radius,
                                                                          height=height)

    return Bz_piecewise_linear_contour


def Bz_circular_single(a, I, cp, calc_radius, height):
    """
    Calculates the Bz field of a single circular coil
    ---------------
    @param a: Turn radius
    @param I: Current in the contour
    @param cp: Calculation domain points
    @param calc_radius: Calculation domain length
    @return: Z-component B of the field of a single circular coil
    """
    mu0 = np.pi * 4e-7
    x = np.linspace(-calc_radius, calc_radius, cp)
    xv, yv = np.meshgrid(x, x)  # Creating meshgrid

    ro = np.sqrt(xv ** 2 + yv ** 2)
    r = np.sqrt(xv ** 2 + yv ** 2 + height ** 2)
    C = mu0 * I / np.pi
    alpha = np.sqrt(a ** 2 + r ** 2 - 2 * a * ro)
    beta = np.sqrt(a ** 2 + r ** 2 + 2 * a * ro)
    k = np.sqrt(1 - alpha ** 2 / beta ** 2)

    Bz = C / (2 * alpha ** 2 * beta) * ((a ** 2 - r ** 2) * ellipe(k ** 2) + alpha ** 2 * ellipk(k ** 2))

    return Bz


def Bz_circular_contour(R, I, spacing, cp, height, a_max):
    """
    Calculates the Bz field for a circular contour
    ---------------
    @param R: Array of radii
    @param I: Current in the contour
    @param spacing: Spacing between coil and the calculation domain boundary
    @param cp: Calculation domain points
    @return: Z-component B of the field of a circular contour
    """
    Bz_circular_contour = np.zeros((cp, cp))

    calc_radius = a_max * spacing

    for r in R:
        Bz_circular_contour += Bz_circular_single(a=r,
                                                  I=I,
                                                  cp=cp,
                                                  calc_radius=calc_radius,
                                                  height=height)
    return Bz_circular_contour


def Bz_square_single(m, n, I, cp, calc_radius, height):
    """
    Calculates the Bz field of a single square coil
    ---------------
    @param m: Side parallel to the x-axis
    @param n: Side parallel to the y-axis
    @param I: Current in the contour
    @param cp: Calculation domain points
    @param calc_radius: Calculation domain length
    @return: Z-component B of the field of single coil
    """
    mu0 = np.pi * 4e-7
    x = np.linspace(-calc_radius, calc_radius, cp)
    xv, yv = np.meshgrid(x, x)  # Creating meshgrid

    C = mu0 * I / (4 * np.pi)

    c1 = xv + m / 2
    c2 = -xv + m / 2
    c3 = -c2
    c4 = -c1

    d1 = yv + n / 2
    d2 = yv + n / 2
    d3 = yv - n / 2
    d4 = yv - n / 2

    r1 = np.sqrt(c1 ** 2 + d1 ** 2 + height ** 2)
    r2 = np.sqrt(c2 ** 2 + d2 ** 2 + height ** 2)
    r3 = np.sqrt(c3 ** 2 + d3 ** 2 + height ** 2)
    r4 = np.sqrt(c4 ** 2 + d4 ** 2 + height ** 2)

    Bz_square = C * (-c2 / (r2 * (r2 + d2)) + d2 / (r2 * (r2 - c2)) - c1 / (r1 * (r1 + d1)) - d1 / (r1 * (r1 + c1)) -
                     c4 / (r4 * (r4 + d4)) + d4 / (r4 * (r4 - c4)) - c3 / (r3 * (r3 + d3)) - d3 / (r3 * (r3 + c3)))

    return Bz_square


def Bz_square_contour(R, X_side, Y_side, I, spacing, cp, height):
    """
    Calculates the Bz field for a square contour
    ---------------
    @param R: Array of radii
    @param X_side: The largest side parallel to the x-axis
    @param Y_side: The largest side parallel to the y-axis
    @param I: Current in the contour
    @param spacing: Spacing between coil and the calculation domain boundary
    @param cp: Calculation domain points
    """
    X_sides, Y_sides = Radii_in_sides_square(R, X_side, Y_side)

    Bz_square_contour = np.zeros((cp, cp))

    calc_radius = 0.5 * max([X_side, Y_side]) * spacing

    for x, y in zip(X_sides, Y_sides):
        Bz_square_contour += Bz_square_single(m=x,
                                              n=y,
                                              I=I,
                                              cp=cp,
                                              calc_radius=calc_radius,
                                              height=height)

    return Bz_square_contour


def eff(Radii_r, Radii_t, d, h):
    """

    @param Radii_r: Receiver radii
    @param Radii_t: Transmitter radii
    @param d: The distance between the centers of the coils
    @param h: The height between the coils
    @return: efficiency coefficient
    """
    M = 0
    for i in range(1, len(Radii_t)+1):
        for j in range(1, len(Radii_r) + 1):
            for n_i in range(1, 361):
                for k_j in range(1, 361):
                    M += np.cos(2 * np.pi / 360 * (n_i - k_j)) /\
                         np.sqrt((d + Radii_t[i] * np.cos(2*np.pi*n_i/360) - Radii_r[j] * np.cos(2*np.pi*k_j/360))**2 +
                                 (Radii_t[i] * np.sin(2*np.pi*n_i/360) - Radii_r[j] * np.sin(2*np.pi*k_j/360))**2 +
                                 h**2)

            M *= Radii_t[i] * Radii_r[j]

    M *= 4 * np.pi**2 * 1e-7 / 360**2
