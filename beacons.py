import numpy as np
from scipy.ndimage.interpolation import shift
from configuration import *
import itertools
from scipy.integrate import simps
from scipy.spatial import KDTree

def bound(value, abs_bound=0.05):
    return max(-abs_bound, min(abs_bound, value))

class Beacon:
    def __init__(self, location, beac_tag):
        self.pt = np.array([location, location])
        self.w = np.array([0., 0.])
        self.v = np.array([[0., 0.],[0., 0.]])

        self.var = [np.NaN, np.NaN]
        self.amp = [np.NaN, np.NaN]
        self.mt = [np.NaN, np.NaN]
        self.ct = [[np.NaN, np.NaN], [np.NaN, np.NaN]]
        self.neigh = []
        self.neigh_toll = []
        # self.neigh_weigh = [[np.NaN], [np.NaN]]
        self.beac_tag = beac_tag

        self.range = clip_range

    def fnc_ants_at_beacon(self, ants):
        # /TODO We don't use this one anymore?
        count = 0
        for ant_tag in ants:
            if self.beac_tag in ants[ant_tag].neigh:
                count += 1
        self.ants_at_beacon = count

    def in_range(self, location):
        return np.linalg.norm(self.pt[1] - location) <= target_range + numeric_step_margin

    def adapt_range(self):
        if adapt_range_option == 'angle':
            if sum(self.v[0]) != 0. and sum(self.v[1]) != 0.:
                angle = np.arccos(np.dot(self.v[0], self.v[1]) / (np.linalg.norm(self.v[0]) *
                                                                  np.linalg.norm(self.v[1]))) * 180 / (np.pi)
                coefficient = (clip_range-min_clip_range)/180
                self.range = coefficient*angle + min_clip_range
        elif adapt_range_option == 'weights':
            if self.w[1] == 0.:
                self.range = clip_range
            elif self.w[0] > 0:
                coefficient = clip_range - min_clip_range
                self.range = -coefficient*min(1,self.w[1]/self.w[0]) + clip_range
            else:
                self.range = min_clip_range

        elif adapt_range_option == 'no adaption':
            self.range = clip_range

    def find_neigh_beacons(self, beacons):
        self.neigh = [beac_tag for beac_tag in beacons if np.linalg.norm(beacons[beac_tag].pt[1]
                                                - self.pt[1]) < beacons[beac_tag].range]
        if not numeric_step_margin:
            self.neigh_toll = self.neigh
        else:
            self.neigh_toll = [beac_tag for beac_tag in beacons.beacons if np.linalg.norm(beacons.beacons[beac_tag].pt[1]
                                                - self.pt[1]) < beacons.beacons[beac_tag].range + numeric_step_margin]

class Beacons:
    def __init__(self, grid, beacon_grid=default_beacon_grid):
        self.grid = grid
        self.beacon_grid = beacon_grid

        self.beacons = dict()

        self.masks = []
        self.map_closest_beacon = []
        self.n = len(self.beacons)

    def non_influenced_points(self):
        influenced_points_dict = dict()
        for beac_tag in self.beacons:
            # influenced_points_dict[beac_tag] = [point for point in self.grid.points if
            #                    (np.linalg.norm(np.array(point) - self.beacons[beac_tag].pt[1]) <= clip_range)]
            influenced_points_dict[beac_tag] = [point for point in self.grid.points if
                (np.linalg.norm(np.array(point) - self.beacons[beac_tag].pt[1]) <= self.beacons[beac_tag].range)]
        influenced_points = list(set([item for l in influenced_points_dict.values() for item in l]))
        non_influenced_points = [point for point in self.grid.points if point not in influenced_points]

        return np.array(non_influenced_points)

    def evaporate_weights(self, rho=default_rho):
        for beac_tag in self.beacons:
            self.beacons[beac_tag].w *= (1 - rho)

    def initialize_weights(self):
        for beac_tag in self.beacons:
            self.beacons[beac_tag].w = np.array([0., 0.])

    def fnc_ants_at_beacons(self, ants):
        for beac_tag in self.beacons:
            self.beacons[beac_tag].fnc_ants_at_beacon(ants)

    def adapt_ranges(self):
        for beac_tag in self.beacons:
            self.beacons[beac_tag].adapt_range()

    @staticmethod
    def extend_mask(mask):
        return ((shift(mask, [1, 1], cval=0) + shift(mask, [1, -1], cval=0) + shift(mask, [-1, 1], cval=0) +
                 shift(mask, [-1, -1], cval=0) + shift(mask, [0, -1], cval=0) + shift(mask, [0, 1], cval=0) +
                 shift(mask, [1, 0], cval=0) + shift(mask, [-1, 0], cval=0)) >= 1) * 1


    def check_weights(self,to_check = 'W1', thres=0.):
        if to_check == 'W1':
            return {beac_tag: self.beacons[beac_tag].w[0] for beac_tag in self.beacons if
                    self.beacons[beac_tag].w[0] > thres}
        elif to_check == 'W2':
            return {beac_tag: self.beacons[beac_tag].w[1] for beac_tag in self.beacons if
                    self.beacons[beac_tag].w[1] > thres}
        elif to_check == 'W':
            return {beac_tag: self.beacons[beac_tag].w[0] + self.beacons[beac_tag].w[1] for beac_tag in self.beacons if
                    self.beacons[beac_tag].w[0] > thres or self.beacons[beac_tag].w[1] > thres}

    def check_ants(self, thres=0):
        return {beac_tag: self.beacons[beac_tag].ants_at_beacon for beac_tag in self.beacons if
                self.beacons[beac_tag].ants_at_beacon > thres}

    def find_neigh_beacons(self):
        for beac_tag in self.beacons:
            self.beacons[beac_tag].find_neigh_beacons(self.beacons)

# test = {beac_tag: simulation.beacons.beacons[beac_tag].w[0] + simulation.beacons.beacons[beac_tag].w[1] for beac_tag in
#         simulation.beacons.beacons if
#         simulation.beacons.beacons[beac_tag].w[0] > 0.00001 or simulation.beacons.beacons[beac_tag].w[1] > 0.00001}
