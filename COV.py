import numpy as np
import math
import copy
from utilities import index_of_element


def if_one(cut, S=0):
    for j in range(len(cut) - 1):
        if cut[j] == 1 and cut[j + 1] == 1:
            S += 1
            return if_one(cut[j + 1: len(cut)], S)
        elif cut[j] == 1 and cut[j + 1] != 1:
            return S


def mask_circular(tiles, r):
    """
    Creates a circular binary mask
    ---------------
    @param tiles: Zero two-dimensional array
    @param r: Radius of the calculation domain is in points
    @return: A circular binary mask
    """
    half_cp = len(tiles) // 2
    for x in range(half_cp - r, half_cp + r):
        for y in range(half_cp - r, half_cp + r):
            if math.sqrt((half_cp - x) ** 2 + (half_cp - y) ** 2) <= r:
                tiles[y][x] = 1


def mask_rectangle(tiles, X_side, Y_side):
    """
    Creates a rectangle binary mask
    ---------------
    @param tiles: Zero two-dimensional array
    @param X_side: Half the side parallel to the x-axis of the calculation domain is in points
    @param Y_side: Half the side parallel to the y-axis of the calculation domain is in points
    @return: A rectangle binary mask
    """
    half_cp = len(tiles) // 2
    for x in range(half_cp - X_side, half_cp + X_side + 1):
        for y in range(half_cp - Y_side, half_cp + Y_side + 1):
            tiles[y][x] = 1


def mask_piecewise_linear(tiles, coords):
    """
    Creates a piecewise linear binary mask
    --------------
    @param tiles: Zero two-dimensional array
    @param coords: Coordinate indexes
    @return: A piecewise linear mask
    """
    for i in range(len(coords)):
        try:
            x1, y1 = coords[i][0], coords[i][1]
            x2, y2 = coords[i + 1][0], coords[i + 1][1]
            X = []
            Y = []
            if y1 > y2:
                for y in range(y2, y1 + 1):
                    Y.append(y)
            elif y1 <= y2:
                for y in range(y1, y2 + 1):
                    Y.append(y)
            if x1 > x2:
                for x in range(x2, x1 + 1):
                    X.append(x)
            elif x1 <= x2:
                for x in range(x1, x2 + 1):
                    X.append(x)
            if len(X) == 1:
                x = x1
                for y in Y:
                    tiles[y][x] = 1
            elif len(Y) == 1:
                y = y1
                for x in X:
                    tiles[y][x] = 1
            else:
                if len(X) > len(Y):
                    Y.clear()
                    k = (y2 - y1) / (x2 - x1)
                    b = y1 - k * x1
                    for x in X:
                        Y.append(round(k * x + b))
                elif len(Y) > len(X):
                    X.clear()
                    k = (y2 - y1) / (x2 - x1)
                    b = y1 - k * x1
                    for y in Y:
                        X.append(round((y - b) / k))
                for x, y in zip(X, Y):
                    tiles[y][x] = 1
            X.clear()
            Y.clear()
        except IndexError:
            x1, y1 = coords[i][0], coords[i][1]
            x2, y2 = coords[0][0], coords[0][1]
            X = []
            Y = []
            if y1 > y2:
                for y in range(y2, y1 + 1):
                    Y.append(y)
            elif y1 <= y2:
                for y in range(y1, y2 + 1):
                    Y.append(y)
            if x1 > x2:
                for x in range(x2, x1 + 1):
                    X.append(x)
            elif x1 <= x2:
                for x in range(x1, x2 + 1):
                    X.append(x)
            if len(X) == 1:
                x = x1
                for y in Y:
                    tiles[y][x] = 1
            elif len(Y) == 1:
                y = y1
                for x in X:
                    tiles[y][x] = 1
            else:
                if len(X) > len(Y):
                    Y.clear()
                    k = (y2 - y1) / (x2 - x1)
                    b = y1 - k * x1
                    for x in X:
                        Y.append(round(k * x + b))
                elif len(Y) > len(X):
                    X.clear()
                    k = (y2 - y1) / (x2 - x1)
                    b = y1 - k * x1
                    for y in Y:
                        X.append(round((y - b) / k))
                for x, y in zip(X, Y):
                    tiles[y][x] = 1
            X.clear()
            Y.clear()

    for y in range(len(tiles)):
        cut = tiles[y, :]
        indexes_contour, indexes_vertexes, indexes_rib = [], [], []
        for x in range(len(cut) - 1):
            if cut[x] == 1:
                indexes_contour.append(x)
                if ([x, y] in coords) and (cut[x + 1] != 1) and (cut[x - 1] != 1):
                    index = index_of_element(coords, [x, y])
                    if (index < len(coords) - 1) and (index > 0):
                        if ((y < coords[index + 1][1]) and (y < coords[index - 1][1])) or (
                                (y > coords[index + 1][1]) and (y > coords[index - 1][1])):
                            indexes_vertexes.append(x)
                    else:
                        if ((y < coords[0][1]) and (y < coords[len(coords) - 1][1])) or (
                                (y > coords[0][1]) and (y > coords[len(coords) - 1][1])):
                            indexes_vertexes.append(x)
                elif ([x, y] in coords) and (cut[x + 1] == 1):
                    last_one = x + 1 + if_one(cut[x + 1:])
                    for i in range(x, last_one):
                        indexes_rib.append(i)
        for element in indexes_vertexes:
            indexes_contour.remove(element)
        for element in indexes_rib:
            indexes_contour.remove(element)

        for index in range(0, len(indexes_contour) - 1, 2):
            for x in range(indexes_contour[index] + 1, indexes_contour[index + 1]):
                tiles[y][x] = 1
        indexes_contour.clear()
        indexes_vertexes.clear()
        indexes_rib.clear()

    for x in range(len(tiles - 1)):
        cut = tiles[:, x]

        for index_y in range(len(cut) - 2):
            if cut[index_y] == 1 and cut[index_y + 1] == 0 and cut[index_y + 2] == 1:
                tiles[index_y + 1][x] = 1
            elif cut[index_y] == 0 and cut[index_y + 1] == 1 and cut[index_y + 2] == 0 and (
                    [x, index_y + 1] not in coords):
                tiles[index_y + 1][x] = 0


def max_min_bz_circular(Bz, max_coil_r, spacing, cp, P):
    calc_radius = max_coil_r * spacing
    cell_size = (2 * calc_radius) / (cp - 1)
    r_cov = round(max_coil_r * np.sqrt(P) / cell_size)  # Uniform area in cells
    half_cp = cp // 2
    res = []
    for x in range(half_cp - r_cov, half_cp + r_cov):
        for y in range(half_cp - r_cov, half_cp + r_cov):
            if math.sqrt((half_cp - x) ** 2 + (half_cp - y) ** 2) <= r_cov:
                res.append(Bz[x][y])
    bz_max = np.max(res)
    bz_min = np.min(res)
    return [bz_min, bz_max]


def max_min_bz_rectangle(Bz, X_side, Y_side, spacing, cp, P):
    calc_radius = 0.5 * max([X_side, Y_side]) * spacing
    cell_size = 2 * calc_radius / (cp - 1)
    X_side_COV = round(X_side * np.sqrt(P) / cell_size) // 2
    Y_side_COV = round(Y_side * np.sqrt(P) / cell_size) // 2
    half_cp = cp // 2
    res = []
    for x in range(half_cp - X_side_COV, half_cp + X_side_COV + 1):
        for y in range(half_cp - Y_side_COV, half_cp + Y_side_COV + 1):
            res.append(Bz[x][y])
    bz_max = np.max(res)
    bz_min = np.min(res)
    return [bz_min, bz_max]


def max_min_bz_piecewise_linear(Bz, coords, spacing, cp, P):
    l = []
    for i in range(len(coords)):
        l.append(np.sqrt((coords[i][0]) ** 2 + (coords[i][1]) ** 2))

    calc_radius = max(l) * spacing
    cell_size = 2 * calc_radius / (cp - 1)
    tiles = np.zeros((cp, cp))

    coords_COV = []
    for point in coords:
        coords_COV.append([round(cp // 2 + point[0] * np.sqrt(P) / cell_size),
                           round(cp // 2 + point[1] * np.sqrt(P) / cell_size)])

    mask_piecewise_linear(tiles, coords_COV)
    Bz_masked = np.multiply(tiles, Bz)
    bz_max = np.max(Bz_masked[np.nonzero(Bz_masked)])
    bz_min = np.min(Bz_masked[np.nonzero(Bz_masked)])
    return [bz_min, bz_max]


def calculation_plane(cell_size, height, cp):
    """"""
    if height < cell_size:
        return int(cp / 2 + 1)
    else:
        return int(height / cell_size + cp / 2)


def COV_circle(Bz, max_coil_r, spacing, P):
    """
    Calculates the coefficient of variation for a circular coil
    --------------
    @param Bz: Field for calculating the COV
    @param max_coil_r: The largest radius of a circular coil [m]
    @param height: Height above the coil [m]
    @param P: The boundary of the calculation of the COV
    @return: COV
    """
    cp = len(Bz)  # Calculation domain

    calc_radius = max_coil_r * spacing  # Calculation domain length
    cell_size = (2 * calc_radius) / (cp - 1)

    tiles = np.zeros((cp, cp))
    r_cov = round(max_coil_r * np.sqrt(P) / cell_size)  # Uniform area in cells
    mask_circular(tiles, r_cov)
    Bz_masked = np.multiply(Bz, tiles)
    B_mean = np.sum(Bz_masked) / np.sum(tiles)
    B_std = np.sqrt(np.sum((Bz_masked - np.multiply(B_mean, tiles)) ** 2) / np.sum(tiles))
    COV = B_std / B_mean
    return COV


def COV_square(Bz, X_side, Y_side, spacing, P):
    """
    Calculates the coefficient of variation for a rectangle coil
    ---------------
    @param Bz: Field for calculating the COV
    @param X_side: Side parallel to the x-axis
    @param Y_side: Side parallel to the y-axis
    @param height: Height above the coil
    @param spacing: Spacing between coil and the calculation domain boundary
    @param P: The boundary of the calculation of the COV
    @return: COV
    """
    cp = len(Bz)

    calc_radius = 0.5 * max([X_side, Y_side]) * spacing
    cell_size = 2 * calc_radius / (cp - 1)
    tiles = np.zeros((cp, cp))

    X_side_COV = round(X_side * np.sqrt(P) / cell_size)
    Y_side_COV = round(Y_side * np.sqrt(P) / cell_size)

    mask_rectangle(tiles, X_side_COV // 2, Y_side_COV // 2)

    Bz_masked = np.multiply(Bz, tiles)

    Bz_mean = np.sum(Bz_masked) / np.sum(tiles)
    Bz_std = np.sqrt(np.sum((np.multiply(Bz_mean, tiles) - Bz_masked) ** 2) / (np.sum(tiles)))

    COV = Bz_std / Bz_mean

    return COV


def COV_piecewise_linear(Bz, coords, spacing, P):
    """
    Calculates the coefficient of variation for a square coil
    ---------------
    @param Bz: Field for calculating the COV
    @param coords: Coordinates of the vertices of a piecewise linear coil [m]
    @param height: Height above the coil [m]
    @param spacing: Spacing between coil and the calculation domain boundary
    @param P: The boundary of the calculation of the COV
    @return: COV
    """
    cp = len(Bz)

    l = []
    for i in range(len(coords)):
        l.append(np.sqrt((coords[i][0]) ** 2 + (coords[i][1]) ** 2))

    calc_radius = max(l) * spacing
    cell_size = 2 * calc_radius / (cp - 1)
    tiles = np.zeros((cp, cp))

    coords_COV = []
    for point in coords:
        coords_COV.append([round(cp // 2 + point[0] * np.sqrt(P) / cell_size),
                           round(cp // 2 + point[1] * np.sqrt(P) / cell_size)])

    mask_piecewise_linear(tiles, coords_COV)
    Bz_masked = np.multiply(tiles, Bz)
    Bz_mean = np.sum(Bz_masked) / np.sum(tiles)
    Bz_std = np.sqrt((np.sum((Bz_masked - np.multiply(Bz_mean, tiles)) ** 2)) / (np.sum(tiles)))
    COV = Bz_std / Bz_mean

    return COV
