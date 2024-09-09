def prop_coeff(R):
    """
    Calculates the radius reduction factor
    ---------------
    """
    R.sort()
    R.reverse()
    prop = []
    for i in range(len(R)):
        prop.append(R[i]/max(R))

    return prop


def Radii_in_sides_square(R, X_side, Y_side, split=False):
    """
    Converts radii to sides of a rectangular coil
    ---------------
    """
    if split:
        split_sides = []
        for coil in R:
            prop = prop_coeff(coil)
            sides = []
            for k in prop:
                sides.append([X_side*k, Y_side*k])
            split_sides.append(sides)
        return split_sides
    else:
        prop = prop_coeff(R=R)
        X_sides, Y_sides = [], []
        for k in prop:
            X_sides.append(X_side * k)
            Y_sides.append(Y_side * k)

        return X_sides, Y_sides


def Radii_in_coords(R, coords_max, split=False):
    """
    Converts radii to vertex coordinates of a piecewise linear coil
    ---------------
    """
    if split:
        split_list_of_coords = []
        for coil in R:
            prop = prop_coeff(R=coil)
            list_of_coords = []
            for k in prop:
                new_coords = []
                for point in coords_max:
                    new_coords.append([point[0] * k, point[1] * k])
                list_of_coords.append(new_coords)
            split_list_of_coords.append(list_of_coords)
        return split_list_of_coords
    else:
        prop = prop_coeff(R=R)
        list_of_coords = []
        for k in prop:
            new_coords = []
            for point in coords_max:
                new_coords.append([point[0] * k, point[1] * k])
            list_of_coords.append(new_coords)
        return list_of_coords


def index_of_element(array, element):
    for index in range(len(array)):
        if array[index] == element:
            return index
