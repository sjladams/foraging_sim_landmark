from domain import *
import random

class Ant:
    def __init__(self, nest_location, food_location, ant_tag, epsilon=default_epsilon, target_range=target_range,
                 init_location = []):
        self.nest_location = np.array(nest_location)
        self.food_location = np.array(food_location)
        self.target_range = target_range

        if any(init_location):
            self.nt = np.array([init_location, init_location])
        else:
            self.nt = np.array([self.nest_location, self.nest_location])  # [nt,nt1]

        self.mode = [0, 0]  # [node_t, node_t1] nt1 has been chosen under mode_t
        self.epsilon = epsilon
        self.trips = 0
        self.move = np.array([[0., 0.], [0.,0.]])

        # self.w = np.array([np.NaN, np.NaN])
        self.ant_tag = ant_tag

        self.neigh = []
        self.neigh_toll = []
        self.neigh_ants = []

        self.obs_avoid_mode = False
        self.exploring_mode = True

    def _search_food(self):
        self.mode[0] = self.mode[1]
        self.mode[1] = 0

    def _search_nest(self):
        self.mode[0] = self.mode[1]
        self.mode[1] = 1

    def _reached_food(self):
        return np.linalg.norm(self.nt[1] - self.food_location) <= self.target_range + numeric_step_margin

    def _reached_nest(self):
        return np.linalg.norm(self.nt[1] - self.nest_location) <= self.target_range + numeric_step_margin

    def _pick_direction(self,beacons,ants):
        if self.mode[1] == 0:
            w_type = 1 # w_type = 1
            w_type_wander = 0
        else:
            w_type = 0 # w_type = 0
            w_type_wander = 1
        self.move[0] = self.move[1]

        # derv = drv_gaussian(self.nt[1][0], self.nt[1][1], beacons.beacons, w_type) + \
        #        np.array([x_drv_elips_gaussian(self.nt[1][0], self.nt[1][1]),
        #                 y_drv_elips_gaussian(self.nt[1][0], self.nt[1][1])])

        # self.find_neigh_beacons(beacons)

        # weights = {beac_tag:beacons.beacons[beac_tag].w[w_type] for beac_tag in self.neigh}
        weights = {beac_tag: beacons.beacons[beac_tag].w[w_type] for beac_tag in self.neigh_toll}

        # weights_wander = {beac_tag: beacons.beacons[beac_tag].w[w_type_wander] for beac_tag in self.neigh}
        max_beac_tag = max(weights, key=weights.get, default = 'nothing_in_range')
        # max_beac_tag_wander = max(weights_wander, key=weights_wander.get, default='nothing_in_range')
        derv = np.zeros(2)
        if max_beac_tag != 'nothing_in_range':
            if beacons.beacons[max_beac_tag].w[w_type] > 0.:
                # location_max_beacon = beacons.beacons[max_beac_tag].pt[1]
                # derv = location_max_beacon - self.nt[1]
                # derv = self.normalize(beacons.beacons[max_beac_tag].v[w_type])

                # derv = beacons.beacons[max_beac_tag].v[w_type]
                derv = sum(np.array([weights[tag] * beacons.beacons[tag].v[w_type] for tag in weights]))

                # if np.isnan(derv[0]) or np.isnan(derv[1]):
                #     print('hier zit ie')

            # # /TODO enable speeding up exploring
            # elif self.neigh_ants:
            #     # # repel_vec = self.normalize(sum(np.array([self.normalize(self.nt[1]-ants[ant_tag].nt[1])
            #     # #                              for ant_tag in self.neigh_ants])))
            #     # repel_vec = sum(np.array([self.normalize(self.nt[1] - ants[ant_tag].nt[1])
            #     #                                          for ant_tag in self.neigh_ants]))
            #     #
            #     # if np.isnan(sum(repel_vec)):
            #     #     derv = self.normalize(np.array([random.uniform(-1, 1), random.uniform(-1, 1)]))
            #     # else:
            #     #     derv = self.normalize(repel_vec + np.array([random.uniform(-1, 1), random.uniform(-1, 1)]))
            #
            #     # repel_vec = self.normalize(sum(np.array([self.normalize(self.nt[1]-ants[ant_tag].nt[1])
            #     #                              for ant_tag in self.neigh_ants])))
            #     repel_vec = sum(np.array([self.normalize(self.nt[1] - beacons.beacons[beac_tag].pt[1])
            #                                              for beac_tag in self.neigh]))
            #
            #     if np.isnan(sum(repel_vec)):
            #         derv = self.normalize(np.array([random.uniform(-1, 1), random.uniform(-1, 1)]))
            #     else:
            #         derv = self.normalize(repel_vec + np.array([random.uniform(-1, 1), random.uniform(-1, 1)]))


        # derv = np.array([x_drv_elips_gaussian(self.nt[1][0], self.nt[1][1]),
        #                  y_drv_elips_gaussian(self.nt[1][0], self.nt[1][1])]) + guided_move

        # if np.linalg.norm(derv) < step_threshold or self.epsilon > random.uniform(0,1):

        if max_beac_tag == 'nothing_in_range':
            derv = np.array([random.uniform(-1, 1), random.uniform(-1, 1)])
            self.exploring_mode = True
        elif beacons.beacons[max_beac_tag].w[w_type] < step_threshold or self.epsilon > random.uniform(0, 1):
            derv = np.array([random.uniform(-1,1),random.uniform(-1,1)])
            self.exploring_mode = True
            # derv = np.random.normal(scale=dt, size=(2))
        elif np.isnan(np.linalg.norm(derv)) or np.linalg.norm(derv)<step_threshold:
            derv = np.array([random.uniform(-1, 1), random.uniform(-1, 1)])
            self.exploring_mode = True
        else:
            self.exploring_mode = False

        # if move_type == 'add':
        #    self.move[1] = self.normalize(self.normalize(derv)*dt + self.move[1])*dt
        # elif move_type == 'der':
        #     self.move[1] = self.normalize(derv)*dt
        # elif move_type == 'add_switch':
        if self.mode[0] != self.mode[1]:
            self.move[1] = -self.move[0]
        else:
            if self.exploring_mode:
                self.move[1] = self.normalize(self.normalize(derv)*dt + self.move[1])*dt
            else:
                self.move[1] = self.normalize(derv) * dt
                # self.move[1] = self.normalize(self.normalize(derv) * dt + self.move[1]) * dt

        return self.move[1]


    def step(self,beacons,ants,grid):
        self.nt[0] = self.nt[1]
        self.nt[1], self.obs_avoid_mode = grid.obstacle_avoidance(self.nt[1], self._pick_direction(beacons,ants))

        if np.isnan(self.nt[1][0]) or np.isnan(self.nt[1][1]):
            print('hoe dan')
        if self.mode[1] == 0 and self._reached_food():
            self.trips += 1
            self._search_nest()
        elif self.mode[1] == 1 and self._reached_nest():
            self.trips += 1
            self._search_food()
        else:
            self.mode[0] = self.mode[1]

    @staticmethod
    def normalize(item):
        return item / np.linalg.norm(item)

    def find_closest_beacon(self, beacons):
        neigh_dist = {beac_tag: np.linalg.norm(beacons.beacons[beac_tag].pt[1] - self.nt[1]) for beac_tag
                      in self.neigh}
        self.cl_beac = min(neigh_dist, key=neigh_dist.get, default=None)

    def find_neigh_beacons(self, beacons):
        self.neigh = [beac_tag for beac_tag in beacons.beacons if np.linalg.norm(beacons.beacons[beac_tag].pt[1]
                                                                - self.nt[1]) < clip_range]
        self.neigh_toll = [beac_tag for beac_tag in beacons.beacons if np.linalg.norm(beacons.beacons[beac_tag].pt[1]
                                        - self.nt[1]) < clip_range + numeric_step_margin]

    def find_neigh_ants(self, ants):
        self.neigh_ants = [ant_tag for ant_tag in ants if (np.linalg.norm(ants[ant_tag].nt[1]
                              - self.nt[1]) < clip_range) and ant_tag != self.ant_tag]

class Ants:
    def __init__(self, nest_location, food_location, epsilon=default_epsilon):
        # self.ants = [Ant(nest_node, food_node, epsilon=default_epsilon) for _ in range(0, N)]
        # self.ants = {ant_tag: Ant(nest_location, food_location, ant_tag, epsilon=epsilon)
        #              for ant_tag in range(1, N+1)}
        self.ants = dict()

        self.nest_location = nest_location
        self.food_location = food_location
        self.epsilon = epsilon

    def steps(self, beacons, grid):
        for ant_tag in self.ants:
            self.ants[ant_tag].step(beacons, self.ants, grid)

    def find_closest_beacon(self, beacons):
        for ant_tag in self.ants:
            self.ants[ant_tag].find_closest_beacon(beacons)



    # def update_weights(self, beacons):
    #     for ant_tag in self.ants:
    #         self.ants[ant_tag].update_weights(beacons)

    def release_ants(self, n, beac_tags):
        next_tag = max(list(self.ants.keys()) + beac_tags,default=-1)+1
        for ant_tag in range(next_tag,next_tag+n):
            self.ants[ant_tag] = Ant(self.nest_location, self.food_location, ant_tag, epsilon=self.epsilon)

    def find_neigh_beacons(self, beacons):
        for ant_tag in self.ants:
            self.ants[ant_tag].find_neigh_beacons(beacons)

    def find_neigh_ants(self):
        for ant_tag in self.ants:
            self.ants[ant_tag].find_neigh_ants(self.ants)

    # def ants_mapper(fnc):
    #     def inner(self, beacons):
    #         for ant_tag in self.ants:
    #             self.ants[ant_tag].fnc
    #     return inner
