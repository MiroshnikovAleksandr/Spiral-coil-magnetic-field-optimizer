import matplotlib.pyplot as plt
import numpy as np
from utilities import Radii_in_coords, Radii_in_sides_square


def plot_2d(Bz, height, a_max, spacing, cp, ax, COV, V):
    """
    Draws a 2D graph of the magnetic field with a cut on the X-axis
    """
    calc_radius = a_max * spacing  # Calculation domain length
    x = np.linspace(-calc_radius, calc_radius, cp)
    ax.plot(x * 1e2, Bz[cp // 2, :] * 1e6, color='#000000',
            label=f'COV = {round(COV * 100, 2)} %' + '\n' + f'V = {round(V * 100, 2)} %')
    plt.xlabel('x [cm]')
    plt.ylabel('Bz [uT]')
    plt.title('Bz Field at {} mm height'.format(height * 1e3))
    ax.grid()
    ax.legend()


def plot_3d(Bz, height, a_max, spacing, cp, ax):
    """
    Draws a 3D graph of the magnetic field
    """
    calc_radius = a_max * spacing  # Calculation domain length
    x = np.linspace(-calc_radius, calc_radius, cp)
    xv, yv = np.meshgrid(x, x)  # Creating meshgrid

    ax.plot_surface(xv * 1e2, yv * 1e2, Bz * 1e6, cmap='inferno')
    ax.set_title('Bz Field at {} mm height'.format(height * 1e3))
    plt.xlabel('$x$, cm')
    plt.ylabel('$y$, cm')


def plot_coil(a_max, spacing, R, ax):
    """
    Draws a circular coil
    ---------------
    @param a_max: Radius of the largest turn
    @param spacing: Spacing between coil and the calculation domain boundary
    @param R: Set of radii
    """
    ax.set_xlim((-a_max * spacing, a_max * spacing))
    ax.set_ylim((-a_max * spacing, a_max * spacing))
    for i, radius in enumerate(R):
        circle = plt.Circle((0, 0), radius, fill=False, color='black')

        ax.add_patch(circle)

        plt.xlabel('x [m]')
        plt.ylabel('y [m]')


def plot_square_coil(m_max, n_max, spacing, R, ax):
    """
    Draws a rectangular coil
    ---------------
    @param m: Side parallel to the x-axis of the largest turn
    @param n: Side parallel to the y-axis of the largest turn
    @param spacing: Spacing between coil and the calculation domain boundary
    @param R: Set of radii
    """
    m_i, n_i = Radii_in_sides_square(R, m_max, n_max)

    max_size = 0.5 * np.max([m_max, n_max]) * spacing

    ax.set_xlim((-max_size, max_size))
    ax.set_ylim((-max_size, max_size))

    for i, size_m in enumerate(m_i):
        for j, size_n in enumerate(n_i):
            if i == j:
                rec = plt.Rectangle((-size_m / 2, -size_n / 2), size_m, size_n, fill=False)
                ax.add_patch(rec)

        plt.xlabel('x [m]')
        plt.ylabel('y [m]')


def plot_piecewise_linear_coil(coords_max, spacing, R, ax):
    """
    Draws a piecewise linear coil
    ---------------
    @param coords_max: Coordinates of the largest turn
    @param spacing: Spacing between coil and the calculation domain boundary
    @param R: Set of radii
    """
    l = []
    for i in range(len(coords_max)-1):
        l.append(np.sqrt(coords_max[i][0]**2 + coords_max[i][1]**2))

    r_max = max(l) * spacing

    ax.set_xlim((-r_max, r_max))
    ax.set_ylim((-r_max, r_max))

    list_of_coords = Radii_in_coords(R, coords_max)

    for coords in list_of_coords:
        for i in range(len(coords)):
            try:
                ax.plot([coords[i][0], coords[i + 1][0]], [coords[i][1], coords[i + 1][1]], color='#000000')
            except IndexError:
                ax.plot([coords[i][0], coords[0][0]], [coords[i][1], coords[0][1]], color='#000000')

    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
