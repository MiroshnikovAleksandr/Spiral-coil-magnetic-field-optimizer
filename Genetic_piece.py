import random
import sys
import math
import time

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tomli

import Bz_Field
import COV
import Resistance
import Plot

from deap import base, algorithms
from deap import creator
from deap import tools
from COV import COV_circle, COV_square, COV_piecewise_linear
from utilities import index_of_element, Radii_in_coords

import warnings

warnings.filterwarnings("error")

toolbox = base.Toolbox()  # create toolbox for genetic algorithm

creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)


class GeneticPiece:
    """
    This class describes the genetic algorithm.
    """

    def __init__(self, params):
        # Genetic parameters
        self.hall_of_fame = None
        self.logbook = None
        self.pop = None
        self.no_of_generations = params['gen']['no_of_generations']
        self.len_of_turn = params['gen']['length_of_turn']
        self.population_size = params['gen']['population_size']
        self.probability_of_mutation = params['gen']['probability_of_mutation']
        self.tournSel_k = params['gen']['tournSel_k']
        self.CXPB = params['gen']['CXPB']
        self.MUTPB = params['gen']['MUTPB']
        self.p_one = params['gen']['p_one']
        self.discreteness = params['gen']['discreteness']
        self.fitn_f = params['gen']['fitn_f']

        # Geometric and electrical parameters
        self.figure = 'Piecewise'
        self.coords = params['geom']['coords']
        self.spacing = params['geom']['spacing']
        self.cp = params['geom']['cp']
        self.minimal_gap = params['geom']['minimal_gap']
        self.height = params['geom']['height']
        self.calculation_area = params['geom']['calculation_area']
        self.material = params['geom']['material']
        self.freq = params['geom']['freq']
        self.I = params['geom']['I']

        self.ind_length = None

    def max_side(self):
        """
        Calculates the coefficient for :var:'minimal_gap' based on the coil geometry.
        @param
        @return: float coefficient
        """
        starting_coords = self.coords
        sides = [(starting_coords[i], starting_coords[i + 1]) for i in range(len(starting_coords) - 1)]
        sides.append((starting_coords[-1], starting_coords[0]))
        sorted_sides = sorted(sides,
                              key=lambda x: np.sqrt((x[0][0] - x[1][0]) ** 2 + (x[0][1] - x[1][1]) ** 2),
                              reverse=True)
        max_side = sorted_sides[0]
        k = (max_side[0][1] - max_side[1][1]) / (max_side[0][0] - max_side[1][0])
        bounds_const = np.sqrt(1 + k ** 2) / (
            abs(max_side[0][1] - k * max_side[0][0]))
        self.ind_length = self.discreteness * int(
            (1 - 0.1) // (bounds_const * self.minimal_gap))  # number of digits encoding a chromosome

    def zero_or_one(self):
        """
        Returns 1 with probability p and 0 with probability 1 - p.
        @param
        @return: 0 or 1
        """
        rand = random.random()
        if rand > self.p_one:
            return 0
        else:
            return 1

    def preparation(self):
        """
        Describes and generates basic genetic algorithm objects:
        an individual is a list of random length consisting of 1s and 0s, that encodes
        a sequence of turn radiuses of a spiral coil;
        the population is a set of individuals.
        @return:
        """
        toolbox.register("ZeroOrOne", self.zero_or_one)
        toolbox.register("individual", tools.initRepeat, creator.Individual,
                         toolbox.ZeroOrOne, self.ind_length)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        self.pop = toolbox.population(n=self.population_size)

    def decode(self, ind: list):
        """
        Decodes the individual chromosome into a list of coeffs of sequential coil turns
        @param ind: chromosome
        @return: list of coeffs

        First, the function iterates over a chromosome and sequentially counts the 1s in the chromosome.
        If it encounters a 0 after a 1 it checks whether the amount of 1s is divisible by :var:'self.discreteness', so
        that the chromosome can be split into an integer number of coil turns. The preceding 1s are replaced with
        0s until the amount of 1s is divisible by :var:'self.discreteness'. In the beginning a 0 is appended to
        :var:'ind' so that when the iteration reaches the end it runs a final check. The 0 is deleted after the
        iteration. The resulting 1s positions are recorded in the :var:'ones_placement' list.

        Secondly, every consecutive :var:'self.discreteness' number of 1s are mapped to a specific wire location, or
        a radius. Each :var:'self.discreteness'-d 1 in :var:'ones_placement' starting with the :var:'self.discreteness'//2-d
        is the center of a wire.
        """
        coeffs = []

        # Append a 0 to :var:'ind'
        ind = ind + [0]
        # Counter of 1s
        cnt = 0
        # List of 1s placement
        ones_placement = []  # list of positions of 1s in ::ind::

        # Iterate over :var:'ind'
        for i in range(self.ind_length):
            if ind[i] == 1:
                cnt += 1
                ones_placement.append(i)
            if ind[i] == 1 and ind[i + 1] == 0 and cnt % self.discreteness != 0:
                j = i
                # Replace preceding 1s with 0s until :var:'cnt' is divisible by :var:'self.discreteness'
                while cnt % self.discreteness != 0:
                    ind[j] = 0
                    cnt -= 1
                    ones_placement.pop()
                    j -= 1

        # Remove the 0 from the end of :var:'ind'
        ind.pop()

        # Map each :var:'self.discreteness' number of 1s to a radius
        for i in range(self.discreteness // 2, len(ones_placement), self.discreteness):
            coeffs.append(0.1 + ones_placement[i] * (1 - 0.1) / self.ind_length)
        # If the :var:'coeffs' array ends up being empty, add a 0 to it
        if not coeffs:
            coeffs.append(0)
        return coeffs

    def determine_Bz(self, ind: list):
        """
        Calculates the Z-component of magnetic inductance Bz at the specified height above the coil encoded
        by :var:'ind' in the specified calculation area.
        @param ind: chromosome
        @return: numpy.ndarray(2, 2) Bz
        """
        R = self.decode(ind)
        return Bz_Field.Bz_piecewise_linear_contour(R=R,
                                                    coords=self.coords,
                                                    I=self.I,
                                                    spacing=self.spacing,
                                                    cp=self.cp,
                                                    height=self.height)

    def determine_COV(self, ind: list):
        """
        Calculates the coefficient (COV) of variation of the Bz of the coil, encoded by :var:'ind'.
        @param ind: chromosome
        @return: tuple (cov, )
        """
        try:
            bz = self.determine_Bz(ind)
        except ZeroDivisionError:
            return 1,
        try:
            cov = COV.COV_piecewise_linear(Bz=bz,
                                           coords=self.coords,
                                           spacing=self.spacing,
                                           P=self.calculation_area)
        # Return 1 if COV is incalculable
        except RuntimeWarning:
            cov = 1
        return cov,

    def max_min_bz_ratio(self, ind: list):
        """
        Calculates (max(Bz) - min(Bz)) / max(Bz) for Bz inside the calculation area.
        @param ind: chromosome
        @return: tuple ((max(Bz) - min(Bz)) / max(Bz), )
        """
        try:
            bz = self.determine_Bz(ind)
        except ZeroDivisionError:
            return 1,
        bz_max_min = COV.max_min_bz_piecewise_linear(Bz=bz,
                                                     coords=self.coords,
                                                     spacing=self.spacing,
                                                     cp=self.cp,
                                                     P=self.calculation_area)
        bz_min = bz_max_min[0]
        bz_max = bz_max_min[1]
        try:
            value = abs((bz_max - bz_min) / bz_max)
        except RuntimeWarning:
            value = 1
        return value,

    def check_feasibility(self, ind: list):
        """
        Checks whether the coil encoded by the :var:'ind' has enough coil turns.
        @param ind: chromosome
        @return: boolean
        """
        return len(self.decode(ind)) == 9

    def execution(self):
        """
        Executes the genetic algorithm with selection, mating and mutation.
        Also stores the best individual from the perspective of its objective function in @hall_of_fame@.
        @return: list, containing the COV of the best individual
        """
        # registering objective function with constraint
        if self.fitn_f == 'max_min':
            toolbox.register("evaluate", self.max_min_bz_ratio)
        elif self.fitn_f == 'COV':
            toolbox.register("evaluate", self.determine_COV)
        # toolbox.decorate("evaluate",
        #                  tools.DeltaPenalty(self.check_feasibility, 1.5))  # constraint on the objective function

        # registering basic processes using built-in functions in DEAP
        toolbox.register("select", tools.selTournament, tournsize=self.tournSel_k)  # selection strategy
        toolbox.register("mate", tools.cxTwoPoint)  # strategy for crossover, this classic two point crossover
        toolbox.register("mutate", tools.mutFlipBit,
                         indpb=self.probability_of_mutation)  # mutation strategy with probability of mutation

        # register statistics
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register('Min', np.min)
        stats.register('Max', np.max)
        stats.register('Avg', np.mean)
        stats.register('Std', np.std)

        self.pop, self.logbook = algorithms.eaSimple(self.pop,
                                                     toolbox,
                                                     cxpb=self.CXPB,
                                                     mutpb=self.MUTPB,
                                                     ngen=self.no_of_generations,
                                                     stats=stats,
                                                     verbose=True)

        self.hall_of_fame = tools.HallOfFame(1)
        self.hall_of_fame.update(self.pop)
        return self.decode(self.hall_of_fame[0])

    def show(self):
        """
        Displays the results (statistics plot, best COV, total length of the best individual).
        @return:
        """
        print(self.hall_of_fame[0])
        print(self.decode(self.hall_of_fame[0]))
        print('1 - Bz_min/Bz_max = ' + str(self.max_min_bz_ratio(self.hall_of_fame[0])))
        print('COV = ' + str(self.determine_COV(self.hall_of_fame[0])))


# with open('parameters.toml', 'rb') as f:
#     parameters = tomli.load(f)
# GA = GeneticPiece(parameters)
# GA.max_side()
# GA.preparation()
# GA.execution()
# GA.show()
# fig, ax = plt.subplots()
# Plot.plot_square_coil(m_max=GA.X_side, n_max=GA.Y_side, spacing=GA.spacing,
#                       R=GA.decode(GA.hall_of_fame[0]), ax=ax)
# Plot.plot_2d(Bz=GA.determine_Bz(GA.hall_of_fame[0]),
#              height=GA.height,
#              a_max=0.25,
#              spacing=GA.spacing,
#              cp=GA.cp,
#              ax=ax,
#              COV=1)
# plt.show()
