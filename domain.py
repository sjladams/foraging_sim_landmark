import numpy as np
from configuration import *
from scipy.spatial import KDTree

class Grid:
    def __init__(self,grid_size=default_grid_size,domain=default_domain, obstacle=None):
        self.grid_size = grid_size
        self.domain = domain
        self.obstacle = obstacle

        self.x = np.linspace(0,domain[0],grid_size[0]+1)
        self.y = np.linspace(0,domain[1],grid_size[1]+1)

        self.X, self.Y = np.meshgrid(self.x,self.y)

        self.points = [tuple([cx, cy]) for cx in self.x for cy in self.y]

        self.ref_lines_domain = np.array([[[0, 0], [0, domain[1]]],
                                     [[0, domain[1]], [domain[0], domain[1]]],
                                     [[domain[0], domain[1]], [domain[0], 0]],
                                     [[domain[0], 0], [0, 0]]])

        if self.obstacle:
            self.ref_lines_obstacle = np.array([[obstacle[0], obstacle[3]],
                                           [obstacle[3], obstacle[2]],
                                           [obstacle[2], obstacle[1]],
                                           [obstacle[1], obstacle[0]]])
            self.ref_lines = np.concatenate((self.ref_lines_domain, self.ref_lines_obstacle), axis=0)
        else:
            self.ref_lines = self.ref_lines_domain

        # self.ref_lines = [np.array([[0,0],[0,domain[1]]]),
        #                  np.array([[0,domain[1]],[domain[0],domain[1]]]),
        #                  np.array([[domain[0],domain[1]],[domain[0],0]]),
        #                  np.array([[domain[0],0],[0,0]])]

        # self.W1 = np.zeros(self.X.shape)
        # self.W2 = np.zeros(self.X.shape)

    # def update_graph_weights(self,beacons):
    #     self.W1 = gaussian(self.X, self.Y, beacons.beacons,0)
    #     self.W2 = gaussian(self.X, self.Y, beacons.beacons, 1)
    #     self.W = self.W1 + self.W2 + np.ones(self.W1.shape)*offset

    @staticmethod
    def normalize(item):
        return item / np.linalg.norm(item)

    @staticmethod
    def perpendicular(a):
        b = np.empty_like(a)
        #     b[0] = -a[1]
        #     b[1] = a[0]
        b[0] = a[1]
        b[1] = -a[0]
        return b

    @staticmethod
    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    # def check_in_domain(self, point):
    #     if (point[0] > 0 and point[0] < self.domain[0] and \
    #             point[1] > 0 and point[1] < self.domain[1]):
    #         return True
    #     else:
    #         return False

    @staticmethod
    def check_direction_vectors(vector_1, vector_2):
        unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
        unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
        dot_product = np.dot(unit_vector_1, unit_vector_2)
        return dot_product > 0

    def check_in_domain(self,point):
        if (point[0] > 0 and point[0] < self.domain[0] and point[1] > 0 and point[1] < self.domain[1]):
            if obstacle:
                if (point[0] > self.obstacle[0][0] and point[0] < self.obstacle[3][0] and \
                        point[1] > self.obstacle[0][1] and point[1] < self.obstacle[1][1]):
                    # point within obstacle area:
                    return False
                else:
                    return True
            else:
                return True
        else:
            return False

    # def line_intersection(self, ant_locations):
    #     for count, ref_line in enumerate(self.ref_lines):
    #         xdiff = np.array([ref_line[0][0] - ref_line[1][0], ant_locations[0][0] - ant_locations[1][0]])
    #         ydiff = np.array([ref_line[0][1] - ref_line[1][1], ant_locations[0][1] - ant_locations[1][1]])
    #
    #         div = self.det(xdiff, ydiff)
    #         if div != 0:
    #             d = (self.det(*ref_line), self.det(*ant_locations))
    #             inter = np.array([self.det(d, xdiff) / div, self.det(d, ydiff) / div])
    #             if np.linalg.norm(inter - ant_locations[0]) < np.linalg.norm(ant_locations[0] - ant_locations[1]):
    #                 return count, inter
    #     return None, None

    def line_intersection(self,ant_locations):
        pot_inter = {}
        move_vec = np.round(np.array([ant_locations[1][0] - ant_locations[0][0],
                                      ant_locations[1][1] - ant_locations[0][1]]), 3)

        for count, ref_line in enumerate(self.ref_lines):
            xdiff = np.array([ref_line[0][0] - ref_line[1][0], ant_locations[0][0] - ant_locations[1][0]])
            ydiff = np.array([ref_line[0][1] - ref_line[1][1], ant_locations[0][1] - ant_locations[1][1]])

            div = self.det(xdiff, ydiff)
            if div != 0:
                d = (self.det(*ref_line), self.det(*ant_locations))
                inter = np.array([self.det(d, xdiff) / div, self.det(d, ydiff) / div])

                inter_vec = np.round(np.array([inter[0] - ant_locations[0][0],
                                               inter[1] - ant_locations[0][1]]), 3)

                if self.check_direction_vectors(move_vec, inter_vec):
                    pot_inter[count] = {'inter': inter,
                                        'inter_dist': np.linalg.norm(inter_vec)}

        if pot_inter:
            key_min = min(pot_inter.keys(), key=(lambda k: pot_inter[k]['inter_dist']))
            return key_min, pot_inter[key_min]['inter']
        else:
            return None, None



    # def obstacle_avoidance(self, start_point, move):
    #     if not self.check_in_domain(start_point + move):
    #         index_line, inter = self.line_intersection(np.array([start_point, start_point + move]))
    #         if isinstance(index_line, int):
    #             ref_line = self.ref_lines[index_line]
    #             ref_vec = np.array([ref_line[1][0] - ref_line[0][0], ref_line[1][1] - ref_line[0][1]])
    #             return inter + self.normalize(self.perpendicular(ref_vec)) * (
    #                         np.linalg.norm(move) - np.linalg.norm(start_point - inter)), True
    #         else:
    #             return start_point + move, False
    #     else:
    #         return start_point + move, False

    def obstacle_avoidance(self,start_point, move):
        no_obs_new_point = start_point + move

        if not self.check_in_domain(start_point + move):
            index_line, inter = self.line_intersection(np.array([start_point, start_point + move]))
            if isinstance(index_line, int):
                ref_line = self.ref_lines[index_line]
                ref_vec = np.array([ref_line[1][0] - ref_line[0][0], ref_line[1][1] - ref_line[0][1]])

                per_vec = self.perpendicular(ref_vec)
                if not per_vec[0]:
                    new_point = np.array([no_obs_new_point[0], 2 * inter[1] - move[1] - start_point[1]])
                else:
                    new_point = np.array([2 * inter[0] - move[0] - start_point[0], no_obs_new_point[1]])

                # Check if mirrored point is within region (we do not calculate double bouncing)
                if self.check_in_domain(new_point):
                    return new_point, True
                else:
                    # if mirrored point is not within region of influence, try to move the opposite direction
                    new_point_alt = start_point - move
                    if self.check_in_domain(new_point_alt):
                        print('mirrored location is not feasible, move opposite direction')
                        return new_point_alt, True

                    # if this also doesn't result in an feasible location, stay were you are
                    else:
                        print('robot is stuck, stays at same location')
                        return start_point, True
            else:
                return start_point + move, False
        else:
            return start_point + move, False


def mapper(fnc):
    def inner(x, y, beacons,w_type):
        value = 0
        # for beacon in beacons:
        #     value += fnc(x, y, beacon,w_type)
        for beac_tag in beacons:
            value += fnc(x, y, beacons[beac_tag],w_type)
        return value
    return inner

def shift_clip(value, beacon, w_type):
    # new_value = value - beacon.w[w_type] * np.exp(-(clip_range**2/(2*beacon.var)))
    if beacon.var[w_type]:
        # new_value = value - ampFactor * np.exp(-(clip_range ** 2 / (2 * beacon.var[w_type])))
        new_value = value - beacon.amp[w_type] * np.exp(-(clip_range ** 2 / (2 * beacon.var[w_type])))
    else:
        new_value = np.zeros(value.shape)
    return np.clip(new_value, 0, np.inf)

@mapper
def gaussian(x,y, beacon,w_type):
    # value = beacon.w[w_type] * np.exp(-((x - beacon.pt[1][0]) ** 2 + (y - beacon.pt[1][1]) ** 2) / (2 * beacon.var))
    if beacon.var[w_type]:
        # value = ampFactor * np.exp(-((x - beacon.pt[1][0]) ** 2 + (y - beacon.pt[1][1]) ** 2) / (2 * beacon.var[w_type]))
        value = beacon.amp[w_type] * np.exp(
            -((x - beacon.pt[1][0]) ** 2 + (y - beacon.pt[1][1]) ** 2) / (2 * beacon.var[w_type]))
    else:
        value = np.zeros(x.shape)
    return shift_clip(value,beacon, w_type)

@mapper
def drv_gaussian(x,y, beacon,w_type):
    # value = beacon.w[w_type] * np.exp(-((x - beacon.pt[1][0]) ** 2 + (y - beacon.pt[1][1]) ** 2) / (2 * beacon.var))
    if beacon.var[w_type]:
        # value = ampFactor * np.exp(
        #     -((x - beacon.pt[1][0]) ** 2 + (y - beacon.pt[1][1]) ** 2) / (2 * beacon.var[w_type]))
        value = beacon.amp[w_type] * np.exp(
            -((x - beacon.pt[1][0]) ** 2 + (y - beacon.pt[1][1]) ** 2) / (2 * beacon.var[w_type]))
        to_return = np.array([-(x - beacon.pt[1][0]) / (beacon.var[w_type]), -(y - beacon.pt[1][1]) / (beacon.var[w_type])]) * shift_clip(
            value, beacon, w_type)
    else:
        to_return = np.zeros(2)
    return to_return


# def elips_gaussian(x,y,a=elips_a,c=elips_c,ampl=elips_ampl):
#     value = ampl*np.exp(-a*(x-default_domain[0]/2)**2 - c*(y-default_domain[1]/2)**2)
#     return np.clip(value,0,0.5)
#
# def x_drv_elips_gaussian(x,y,a=elips_a,c=elips_c,ampl=elips_ampl):
#     value = elips_gaussian(x,y,a=a,c=c,ampl=ampl)
#     if value == 0.5:
#         return 0
#     else:
#         return -2*a*(x-default_domain[0]/2)*value
#
# def y_drv_elips_gaussian(x,y,a=elips_a,c=elips_c,ampl=elips_ampl):
#     value = elips_gaussian(x,y,a=a,c=c,ampl=ampl)
#     if value == 0.5:
#         return 0
#     else:
#         return -2*c*(y-default_domain[1]/2)*value

