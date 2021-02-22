from beacons import *
from configuration import *
from domain import *
from ants import *

import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import simps
from scipy.spatial import KDTree
from scipy.ndimage.interpolation import shift
from scipy.spatial import Voronoi, voronoi_plot_2d

if local:
    FOLDER_LOCATION = './figures/'
else:
    FOLDER_LOCATION = './figures_manuels_desk/'

class Simulations:
    def __init__(self,grid_size=default_grid_size, beacon_grid=default_beacon_grid,
                 nest_location=default_nest_location, food_location=default_food_location, N_total=default_N_total,
                 N_batch = default_N_batch, rho=default_rho,rho_v=default_rho_v, domain=default_domain):
        self.nest_location = default_nest_location
        self.food_location = default_food_location
        self.N_total = N_total
        self.N_batch = N_batch
        self.rho = rho
        self.rho_v = rho_v

        self.grid = Grid(grid_size=grid_size,domain=domain)
        self.total_trips = dict()

        self.beacons = Beacons(self.grid, beacon_grid=beacon_grid)
        # self.beacons.initialize_weights()
        # self.beacons.update_masks()
        # self.beacons.update_neighbours_beacons()

        self.ants = Ants(nest_location, food_location, epsilon=default_epsilon)
        self.ants.release_ants(1,list(self.beacons.beacons.keys()))
        # self.ants.steps(self.beacons, self.grid)

        self.update_ant_beacon_connection()
        self.switch_step()
        self.update_ant_beacon_connection()

        # self.update_beacon_weights()
        # self.beacons.update_neighbours_beacons()

        # self.grid.update_graph_weights(self.beacons)
        # self.beacons.update_m_c_beacons(self.grid.W)

        self.food_found = False

    def sim_step(self, time_step, switch_time=250):
        N_till_now = (time_step+1)*self.N_batch

        # ACTION
        if N_till_now < self.N_total:
            self.ants.release_ants(self.N_batch, list(self.beacons.beacons.keys()))
        # UPDATE
        self.update_ant_beacon_connection()

        # ACTION
        self.ants.steps(self.beacons, self.grid)

        # UPDATE
        self.update_ant_beacon_connection()

        # ACTION
        if time_step >= switch_time:
            self.switch_step()

        # UPDATE
        ## Done within switch_step

        # ACTION
        self.beacons.evaporate_weights(rho=self.rho)
        self.update_beacon_weights()

        # UPDATE
        self.update_ant_beacon_connection()
        # self.beacons.update_beacon_configuration(position_changed=False)

        # self.grid.update_graph_weights(self.beacons)

        # STORE
        self.store_nr_trips(time_step)

    # def sim_step_mov_beac(self,time_step,switch_time=250):
    #     self.sim_step(time_step)
    #
    #     if time_step >= switch_time:
    #         # self.beacons.move_step(self.grid.W)
    #         self.switch_step()
    #         # self.grid.update_graph_weights(self.beacons)
    #
    #     self.store_nr_trips(time_step)

    def update_ant_beacon_connection(self):
        self.ants.find_neigh_beacons(self.beacons)           # Depends on location of beacons and ants
        self.ants.find_neigh_ants()  # Depends on location of beacons and ants
        self.ants.find_closest_beacon(self.beacons)                 # Depends on location of beacons and ants, and neigh_beacons
        self.beacons.fnc_ants_at_beacons(self.ants.ants)            # Depends on cl_beaon

    def switch_step(self):
        # self.ants.update_weights(self.beacons)

        old_ants = self.ants.ants.copy()
        tags_changed = []
        for ant_tag in old_ants:
            # if sum(self.ants.ants[ant_tag].w) <= 0.:
            if self.ants.ants[ant_tag].cl_beac == None or self.ants.ants[ant_tag]._reached_food() and not self.food_found:
                # \TODO GEBEUND, AANPASSEN?
                if self.ants.ants[ant_tag]._reached_food() and not self.food_found:
                    self.food_found = True
                    self.food_location = self.ants.ants[ant_tag].nt[1]
            # if self.ants.ants[ant_tag].cl_beac == None:

                tags_changed += [ant_tag]

                self.beacons.beacons[ant_tag] = Beacon(self.ants.ants[ant_tag].nt[1], ant_tag)

                self.initialize_beacon_weights(ant_tag) # only one agents initializes, so no need to update nr of ants per beacon
                # self.beacons.update_beacon_configuration(position_changed=False)

                del self.ants.ants[ant_tag]
                self.update_ant_beacon_connection()

                # # check beacons one by one:
                # self.ants.update_weights(self.beacons)

        weight_dict = self.beacons.check_weights(to_check = 'W',thres=threshold)
        ant_dict = self.beacons.check_ants()

        # /TODO switch back
        old_beacons = self.beacons.beacons.copy()
        for beac_tag in old_beacons:
            if beac_tag not in weight_dict and beac_tag not in ant_dict and beac_tag not in tags_changed:
                self.ants.ants[beac_tag] = Ant(self.nest_location, self.food_location, beac_tag,
                                               init_location=self.beacons.beacons[beac_tag].pt[1])
                del self.beacons.beacons[beac_tag]

                self.update_ant_beacon_connection()


    def reward(self, weights,rew):
        # return self.rho * (lam * max(weights, default=0) + rew) / (ants_at_beacon)
        return self.rho * (lam * max(weights, default=0) + rew)

    def initialize_beacon_weights(self,tag):
        # /TODO remove initialization

        self.ants.ants[tag].find_neigh_beacons(self.beacons)
        W1_weights = [self.beacons.beacons[beac_tag].w[0] for beac_tag in self.ants.ants[tag].neigh]
        W2_weights = [self.beacons.beacons[beac_tag].w[1] for beac_tag in self.ants.ants[tag].neigh]

        ant_w_update = np.array([0., 0.])
        beac_v = np.array([[0., 0.], [0., 0.]])

        if self.ants.ants[tag].mode[0] == 0:
            if self.ants.ants[tag]._reached_nest():
                ant_w_update[0] += self.reward(W1_weights, rew)

                if use_weights_updating_v:
                    beac_v[0] += - self.ants.ants[tag].move[1] * ant_w_update[0]
                else:
                    beac_v[0] += - self.ants.ants[tag].move[1]

            elif self.ants.ants[tag]._reached_food():
                ant_w_update[1] += self.reward(W2_weights, rew)
                if use_weights_updating_v:
                    beac_v[1] += self.ants.ants[tag].move[1] * ant_w_update[1]
                else:
                    beac_v[1] += self.ants.ants[tag].move[1]

            else:
                ant_w_update[0] += self.reward(W1_weights, 0)
                # / TODO check this trange condition
                if self.ants.ants[tag].mode[0] == self.ants.ants[tag].mode[1]:
                    if use_weights_updating_v:
                        beac_v[0] += -self.ants.ants[tag].move[1] * ant_w_update[0]
                    else:
                        beac_v[0] += -self.ants.ants[tag].move[1]


        elif self.ants.ants[tag].mode[0] == 1:
            if self.ants.ants[tag]._reached_food():
                ant_w_update[1] += self.reward(W2_weights, rew)

                if use_weights_updating_v:
                    beac_v[1] += - self.ants.ants[tag].move[1] * ant_w_update[1]
                else:
                    beac_v[1] += - self.ants.ants[tag].move[1]

            elif self.ants.ants[tag]._reached_nest():
                ant_w_update[0] += self.reward(W1_weights, rew)
                if use_weights_updating_v:
                    beac_v[0] += self.ants.ants[tag].move[1] * ant_w_update[0]
                else:
                    beac_v[0] += self.ants.ants[tag].move[1]

            else:
                ant_w_update[1] += self.reward(W2_weights, 0)
                if self.ants.ants[tag].mode[0] == self.ants.ants[tag].mode[1]:
                    if use_weights_updating_v:
                        beac_v[1] += - self.ants.ants[tag].move[1] * ant_w_update[1]
                    else:
                        beac_v[1] += - self.ants.ants[tag].move[1]


        self.beacons.beacons[tag].w[0] += ant_w_update[0]
        self.beacons.beacons[tag].w[1] += ant_w_update[1]

        if use_weights_updating_v:
            if use_rhov_2_init:
                self.beacons.beacons[tag].v[0] += self.rho_v * beac_v[0] / ant_w_update[0]
                self.beacons.beacons[tag].v[1] += self.rho_v * beac_v[1] / ant_w_update[1]
            else:
                self.beacons.beacons[tag].v[0] += beac_v[0] / ant_w_update[0]
                self.beacons.beacons[tag].v[1] += beac_v[1] / ant_w_update[1]
        else:
            if use_rhov_2_init:
                self.beacons.beacons[tag].v[0] += self.rho_v * beac_v[0]
                self.beacons.beacons[tag].v[1] += self.rho_v * beac_v[1]
            else:
                self.beacons.beacons[tag].v[0] += beac_v[0]
                self.beacons.beacons[tag].v[1] += beac_v[1]

        # if self.beacons.beacons[tag].nest_in_range(self.nest_location):
        #     self.beacons.beacons[tag].v[0] = np.array([0,0])
        # elif self.beacons.beacons[tag].food_in_range(self.food_location):
        #     self.beacons.beacons[tag].v[1] = np.array([0, 0])

        if (self.beacons.beacons[tag].pt[1] == self.nest_location).all():
            self.beacons.beacons[tag].v[0] = np.array([0,0])
        elif (self.beacons.beacons[tag].pt[1] == self.food_location).all():
            self.beacons.beacons[tag].v[1] = np.array([0, 0])

        # if (np.linalg.norm(self.beacons.beacons[tag].v[1]) >dt + 0.001 or \
        #         np.linalg.norm(self.beacons.beacons[tag].v[1]) < dt - 0.001) and np.linalg.norm(self.beacons.beacons[tag].v[1]) != 0.:
        #     print('hoe kan dit')
        # if (np.linalg.norm(self.beacons.beacons[tag].v[0]) >dt + 0.001 or \
        #         np.linalg.norm(self.beacons.beacons[tag].v[0]) < dt - 0.001) and np.linalg.norm(self.beacons.beacons[tag].v[0]) != 0.:
        #     print('hoe kan dit')

    def update_beacon_weights(self):
        # /TODO We now reward for every ant that finds the food! Not in line with our concept
        # /TODO check if it is correct to use mode[0]

        for beac_tag in self.beacons.beacons:
            self.beacons.beacons[beac_tag].v *= (1 - self.rho_v)

            ants_at_beacon = [ant_tag for ant_tag in self.ants.ants if beac_tag in self.ants.ants[ant_tag].neigh_toll]

            # if self.beacons.beacons[beac_tag].ants_at_beacon == 0.:
            if not len(ants_at_beacon):
                continue

            beac_v = np.array([[0.,0.],[0.,0.]])
            count_v0 = 0
            count_v1 = 0
            ants_w_update = np.array([0.,0.])

            for ant_tag in ants_at_beacon:
                if self.ants.ants[ant_tag].obs_avoid_mode:
                    continue

                W1_weights = [self.beacons.beacons[beac_tag].w[0] for beac_tag in self.ants.ants[ant_tag].neigh]
                W2_weights = [self.beacons.beacons[beac_tag].w[1] for beac_tag in self.ants.ants[ant_tag].neigh]
                ant_w_update = np.array([0.,0.])

                if self.ants.ants[ant_tag].mode[0] == 0:
                    if self.ants.ants[ant_tag]._reached_nest():
                        ant_w_update[0] = self.reward(W1_weights,rew)

                        count_v0 += 1
                        if use_weights_updating_v:
                            beac_v[0] += -self.ants.ants[ant_tag].move[1] * ant_w_update[0]
                        else:
                            beac_v[0] += -self.ants.ants[ant_tag].move[1]

                    elif self.ants.ants[ant_tag]._reached_food():
                        ant_w_update[1] = self.reward(W2_weights, rew)

                        count_v1 += 1
                        if use_weights_updating_v:
                            beac_v[1] += self.ants.ants[ant_tag].move[1] * ant_w_update[1]
                        else:
                            beac_v[1] += self.ants.ants[ant_tag].move[1]
                    else:
                        ant_w_update[0] = self.reward(W1_weights,0,)

                        count_v0 += 1
                        if use_weights_updating_v:
                            beac_v[0] += -self.ants.ants[ant_tag].move[1] * ant_w_update[0]
                        else:
                            beac_v[0] += -self.ants.ants[ant_tag].move[1]

                elif self.ants.ants[ant_tag].mode[0] == 1:
                    if self.ants.ants[ant_tag]._reached_food():
                        ant_w_update[1] = self.reward(W2_weights,rew)

                        count_v1 += 1
                        if use_weights_updating_v:
                            beac_v[1] += -self.ants.ants[ant_tag].move[1]*ant_w_update[1]
                        else:
                            beac_v[1] += -self.ants.ants[ant_tag].move[1]

                    elif self.ants.ants[ant_tag]._reached_nest():
                        ant_w_update[0] = self.reward(W1_weights,rew)

                        count_v0 += 1
                        if use_weights_updating_v:
                            beac_v[0] += self.ants.ants[ant_tag].move[1]*ant_w_update[0]
                        else:
                            beac_v[0] += self.ants.ants[ant_tag].move[1]
                    else:
                        ant_w_update[1] = self.reward(W2_weights, 0)

                        count_v1 +=1
                        if use_weights_updating_v:
                            beac_v[1] += -self.ants.ants[ant_tag].move[1]*ant_w_update[1]
                        else:
                            beac_v[1] += -self.ants.ants[ant_tag].move[1]

                ants_w_update += ant_w_update
                # self.beacons.beacons[beac_tag].w[0] += ant_w_update[0] / self.beacons.beacons[beac_tag].ants_at_beacon
                self.beacons.beacons[beac_tag].w[0] += ant_w_update[0] / len(ants_at_beacon)
                # self.beacons.beacons[beac_tag].w[1] += ant_w_update[1] /self.beacons.beacons[beac_tag].ants_at_beacon
                self.beacons.beacons[beac_tag].w[1] += ant_w_update[1] / len(ants_at_beacon)


            if np.linalg.norm(self.beacons.beacons[beac_tag].v[1]) and count_v1:
                if use_weights_updating_v:
                    self.beacons.beacons[beac_tag].v[1] += self.rho_v * beac_v[1] / (count_v1* ants_w_update[1])
                else:
                    self.beacons.beacons[beac_tag].v[1] += self.rho_v * beac_v[1] / (count_v1)
            elif count_v1:
                if use_weights_updating_v:
                    if use_rhov_2_init:
                        self.beacons.beacons[beac_tag].v[1] += self.rho_v * beac_v[1] /(count_v1* ants_w_update[1])
                    else:
                        self.beacons.beacons[beac_tag].v[1] += beac_v[1] / (count_v1* ants_w_update[1])
                else:
                    if use_rhov_2_init:
                        self.beacons.beacons[beac_tag].v[1] += self.rho_v * beac_v[1] / (count_v1)
                    else:
                        self.beacons.beacons[beac_tag].v[1] += beac_v[1] / (count_v1)

            if np.linalg.norm(self.beacons.beacons[beac_tag].v[0]) and count_v0:
                if use_weights_updating_v:
                    self.beacons.beacons[beac_tag].v[0] += self.rho_v * beac_v[0] / (count_v0* ants_w_update[0])
                else:
                    self.beacons.beacons[beac_tag].v[0] += self.rho_v * beac_v[0] / (count_v0)
            elif count_v0:
                if use_weights_updating_v:
                    if use_rhov_2_init:
                        self.beacons.beacons[beac_tag].v[0] += self.rho_v * beac_v[0] /(count_v0* ants_w_update[0])
                    else:
                        self.beacons.beacons[beac_tag].v[0] += beac_v[0] / (count_v0* ants_w_update[0])
                else:
                    if use_rhov_2_init:
                        self.beacons.beacons[beac_tag].v[0] += self.rho_v * beac_v[0] / (count_v0)
                    else:
                        self.beacons.beacons[beac_tag].v[0] += beac_v[0] / (count_v0)

            if self.beacons.beacons[beac_tag].nest_in_range(self.nest_location):
                self.beacons.beacons[beac_tag].v[0] = np.array([0, 0])
            elif self.beacons.beacons[beac_tag].food_in_range(self.food_location):
                self.beacons.beacons[beac_tag].v[1] = np.array([0, 0])

            # if (self.beacons.beacons[beac_tag].pt[1] == self.nest_location).all():
            #     self.beacons.beacons[beac_tag].v[0] = np.array([0, 0])
            # elif (self.beacons.beacons[beac_tag].pt[1] == self.food_location).all():
            #     self.beacons.beacons[beac_tag].v[1] = np.array([0, 0])

            # if (np.linalg.norm(self.beacons.beacons[beac_tag].v[1]) > dt + 0.01 or \
            #         np.linalg.norm(self.beacons.beacons[beac_tag].v[1]) < dt - 0.01) and np.linalg.norm(self.beacons.beacons[beac_tag].v[1]) != 0.:
            #     print('hoe kan dit')
            # if (np.linalg.norm(self.beacons.beacons[beac_tag].v[0]) > dt + 0.01 or \
            #         np.linalg.norm(self.beacons.beacons[beac_tag].v[0]) < dt - 0.01) and np.linalg.norm(self.beacons.beacons[beac_tag].v[0]) != 0.:
            #     print('hoe kan dit')

        #
        # for ant_tag in self.ants.ants:
        #     if self.ants.ants[ant_tag].cl_beac == None:
        #         test =1
        #         continue
        #     # self.ants.ants[ant_tag].find_neigh_beacons(self.beacons)
        #     W1_weights = [self.beacons.beacons[beac_tag].w[0] for beac_tag in self.ants.ants[ant_tag].neigh]
        #     # W1_weights = self.ants.ants[ant_tag].neigh_weigh[0].values()
        #     W2_weights = [self.beacons.beacons[beac_tag].w[1] for beac_tag in self.ants.ants[ant_tag].neigh]
        #     # W2_weights = self.ants.ants[ant_tag].neigh_weigh[1].values()
        #
        #     if self.ants.ants[ant_tag].mode[0]==0:
        #         if self.ants.ants[ant_tag]._reached_nest():
        #             self.beacons.beacons[self.ants.ants[ant_tag].cl_beac].w[0] += self.reward(W1_weights,
        #                                                 rew, self.beacons.beacons[self.ants.ants[ant_tag].cl_beac].ants_at_beacon)
        #         elif self.ants.ants[ant_tag]._reached_food():
        #             self.beacons.beacons[self.ants.ants[ant_tag].cl_beac].w[1] += self.reward(W2_weights,
        #                                                 rew, self.beacons.beacons[self.ants.ants[ant_tag].cl_beac].ants_at_beacon)
        #         else:
        #             self.beacons.beacons[self.ants.ants[ant_tag].cl_beac].w[0] += self.reward(W1_weights,
        #                                                 0, self.beacons.beacons[self.ants.ants[ant_tag].cl_beac].ants_at_beacon)
        #     elif self.ants.ants[ant_tag].mode[0]==1:
        #         if self.ants.ants[ant_tag]._reached_food():
        #             self.beacons.beacons[self.ants.ants[ant_tag].cl_beac].w[1] += self.reward(W2_weights,
        #                                                 rew, self.beacons.beacons[self.ants.ants[ant_tag].cl_beac].ants_at_beacon)
        #         elif self.ants.ants[ant_tag]._reached_nest():
        #             self.beacons.beacons[self.ants.ants[ant_tag].cl_beac].w[0] += self.reward(W1_weights,
        #                                                 rew, self.beacons.beacons[self.ants.ants[ant_tag].cl_beac].ants_at_beacon)
        #         else:
        #             self.beacons.beacons[self.ants.ants[ant_tag].cl_beac].w[1] += self.reward(W2_weights,
        #                                                 0, self.beacons.beacons[self.ants.ants[ant_tag].cl_beac].ants_at_beacon)


    def plt(self, to_plot='W'):
        # vor = Voronoi([item.pt[1] for item in self.beacons.beacons])
        # voronoi_plot_2d(vor, show_vertices=False)

        if to_plot == 'W1':
            plt.contourf(self.grid.X, self.grid.Y, self.grid.W1)  # ,levels=np.linspace(0,1,10))
        elif to_plot == 'W2':
            plt.contourf(self.grid.X, self.grid.Y, self.grid.W2)  # ,levels=np.linspace(0,1,10))
        elif to_plot == 'W':
            plt.contourf(self.grid.X, self.grid.Y, self.grid.W)  # ,levels=np.linspace(0,1,10))

        # plt.plot([item.pt[1][0] for item in self.beacons.beacons],
        #          [item.pt[1][1] for item in self.beacons.beacons], 'r*')
        plt.plot([item.nt[1][0] for item in self.ants.ants if item.mode[1] == 0],
                 [item.nt[1][1] for item in self.ants.ants if item.mode[1] == 0], 'g*')
        plt.plot([item.nt[1][0] for item in self.ants.ants if item.mode[1] == 1],
                 [item.nt[1][1] for item in self.ants.ants if item.mode[1] == 1], 'y*')
        plt.plot([self.nest_location[0], self.food_location[0]],
                 [self.nest_location[1], self.food_location[1]], 'r*')
        # plt.plot(list(itertools.chain.from_iterable(self.grid.X)),
        #          list(itertools.chain.from_iterable(self.grid.Y)), 'b*')

        plt.xlim(0-5, self.grid.domain[0]+5)
        plt.ylim(0-5, self.grid.domain[1]+5)
        plt.colorbar()
        plt.show()

    def plt_3d(self, to_plot='W', fig_tag=None):
        fig = plt.figure(figsize=(12, 6))
        ax = fig.gca(projection='3d')

        locations_non_influenced_points = self.beacons.non_influenced_points()
        value_non_influenced_points = np.zeros(locations_non_influenced_points.shape[0])

        locations_beacons = np.array([self.beacons.beacons[beac_tag].pt[1] for beac_tag in self.beacons.beacons])

        if to_plot == 'W1':
            w_beacons = np.array([self.beacons.beacons[beac_tag].w[0] for beac_tag in self.beacons.beacons])
        elif to_plot == 'W2':
            w_beacons = np.array([self.beacons.beacons[beac_tag].w[1] for beac_tag in self.beacons.beacons])
        else:
            w_beacons = np.array([self.beacons.beacons[beac_tag].w[0] + self.beacons.beacons[beac_tag].w[1]
                         for beac_tag in self.beacons.beacons])

        W_beacons = griddata(np.concatenate((locations_beacons,locations_non_influenced_points)),
                             np.concatenate((w_beacons,value_non_influenced_points)),  (self.grid.X, self.grid.Y), method='linear')

        ax.plot_surface(self.grid.X, self.grid.Y, W_beacons,
                        cmap=cm.coolwarm,
                        linewidth=0,
                        antialiased=True)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        if to_plot == 'W1' and fig_tag:
            plt.savefig(FOLDER_LOCATION + 'W1_3d/' + str(to_plot) + '_' + str(fig_tag) + '.png')
            plt.close()
        elif to_plot == 'W2' and fig_tag:
            plt.savefig(FOLDER_LOCATION + 'W2_3d/' + str(to_plot) + '_' + str(fig_tag) + '.png')
            plt.close()
        elif to_plot == 'W' and fig_tag:
            plt.savefig(FOLDER_LOCATION + 'W_3d/' + str(to_plot) + '_' + str(fig_tag) + '.png')
            plt.close()
        else:
            plt.show()
            plt.close()

    def plt_beacons(self, to_plot='W', fig_tag=None):
        # vor = Voronoi([item.pt[1] for item in self.beacons.beacons])

        for beac_tag in self.beacons.beacons:
            if np.isnan(self.beacons.beacons[beac_tag].pt[1][0]) or np.isnan(self.beacons.beacons[beac_tag].pt[1][1]):
                print('hoe kan dit nou weer?')

        if len(self.beacons.beacons) >3:
            vor = Voronoi([self.beacons.beacons[beac_tag].pt[1] for beac_tag in self.beacons.beacons])
            voronoi_plot_2d(vor, show_vertices=False)

        if to_plot == 'W1':
            for beac_tag in self.beacons.check_weights(to_check='W1'):
                item = self.beacons.beacons[beac_tag]
                size = (item.w[0] / max(self.beacons.check_weights(to_check='W1').values())) * 10
                plt.plot([item.pt[1][0]], [item.pt[1][1]], 'o', color='black', markersize=size)
        elif to_plot == 'W2':
            for beac_tag in self.beacons.check_weights(to_check='W2'):
                item = self.beacons.beacons[beac_tag]
                size = (item.w[1] / max(self.beacons.check_weights(to_check='W2').values())) * 10
                plt.plot([item.pt[1][0]], [item.pt[1][1]], 'o', color='black', markersize=size)
        elif to_plot == 'W':
            for beac_tag in self.beacons.check_weights(to_check='W'):
                item = self.beacons.beacons[beac_tag]
                size = ((item.w[1] +item.w[0] ) / max(self.beacons.check_weights(to_check='W').values())) * 10
                plt.plot([item.pt[1][0]], [item.pt[1][1]], 'o', color='black', markersize=size)


        plt.plot([self.nest_location[0], self.food_location[0]],
                 [self.nest_location[1], self.food_location[1]], 'r*')

        plt.plot([self.ants.ants[ant_tag].nt[1][0] for ant_tag in self.ants.ants if
                  self.ants.ants[ant_tag].mode[1] == 0],
                 [self.ants.ants[ant_tag].nt[1][1] for ant_tag in self.ants.ants if
                  self.ants.ants[ant_tag].mode[1] == 0], 'g*')
        plt.plot([self.ants.ants[ant_tag].nt[1][0] for ant_tag in self.ants.ants if
                  self.ants.ants[ant_tag].mode[1] == 1],
                 [self.ants.ants[ant_tag].nt[1][1] for ant_tag in self.ants.ants if
                  self.ants.ants[ant_tag].mode[1] == 1], 'y*')

        # plt.xlim(0, self.grid.domain[0])
        plt.xlim(-2, self.grid.domain[0]+2)
        # plt.ylim(0, self.grid.domain[1])
        plt.ylim(-2, self.grid.domain[1]+2)
        # plt.colorbar()
        if to_plot == 'W1' and fig_tag:
            plt.savefig(FOLDER_LOCATION + 'W1/' + str(to_plot) + '_' + str(fig_tag) + '.png')
            plt.close()
        elif to_plot == 'W2' and fig_tag:
            plt.savefig(FOLDER_LOCATION +'W2/' + str(to_plot) + '_' + str(fig_tag) + '.png')
            plt.close()
        elif to_plot == 'W' and fig_tag:
            plt.savefig(FOLDER_LOCATION +'W/' + str(to_plot) + '_' + str(fig_tag) + '.png')
            plt.close()
        else:
            plt.show()
            plt.close()


    def plt_only_weights_and_vectors(self, to_plot='W', fig_tag=None):
        if len(self.beacons.beacons) >3:
            vor = Voronoi([self.beacons.beacons[beac_tag].pt[1] for beac_tag in self.beacons.beacons])
            voronoi_plot_2d(vor, show_vertices=False)

        if to_plot == 'W1':
            # /TODO change this to all beacons
            for beac_tag in self.beacons.check_weights(to_check='W1'):
                item = self.beacons.beacons[beac_tag]
                # size = (item.w[0] / max(self.beacons.check_weights(to_check='W1').values())) * 10
                size = np.log(item.w[0] / max(self.beacons.check_weights(to_check='W1').values()) + 1) * 10

                plt.plot([item.pt[1][0]], [item.pt[1][1]], 'o', color='black', markersize=size)

                if np.linalg.norm(item.v[0]) >step_threshold and item.w[0] > step_threshold:
                    arrow = self.normalize(item.v[0]) * dt
                    plt.plot([item.pt[1][0],item.pt[1][0]+arrow[0]], [item.pt[1][1],item.pt[1][1]+arrow[1]], color='black')

                # if np.linalg.norm(np.array([item.v[0][0],item.v[0][1]])) > 0.55 or \
                #         np.linalg.norm(np.array([item.v[0][0],item.v[0][1]])) < 0.45:
                #     print('aaap')

        elif to_plot == 'W2':
            for beac_tag in self.beacons.check_weights(to_check='W2'):
                item = self.beacons.beacons[beac_tag]
                # size = (item.w[1] / max(self.beacons.check_weights(to_check='W2').values())) * 10
                size = np.log(item.w[1] / max(self.beacons.check_weights(to_check='W2').values()) + 1) * 10
                plt.plot([item.pt[1][0]], [item.pt[1][1]], 'o', color='black', markersize=size)

                if np.linalg.norm(item.v[1]) >step_threshold and item.w[1] > step_threshold:
                    arrow = self.normalize(item.v[1]) * dt
                    plt.plot([item.pt[1][0],item.pt[1][0]+arrow[0]], [item.pt[1][1],item.pt[1][1]+arrow[1]], color='black')

        elif to_plot == 'W':
            for beac_tag in self.beacons.check_weights(to_check='W'):
                item = self.beacons.beacons[beac_tag]
                # size = ((item.w[1] +item.w[0] ) / max(self.beacons.check_weights(to_check='W').values())) * 10
                size = np.log((item.w[1] + item.w[0]) / max(self.beacons.check_weights(to_check='W').values())+1) * 10
                plt.plot([item.pt[1][0]], [item.pt[1][1]], 'o', color='black', markersize=size)

                if np.linalg.norm(item.v[0]) >step_threshold and item.w[0] > step_threshold:
                    arrow0 = self.normalize(item.v[0]) *dt
                    plt.plot([item.pt[1][0], item.pt[1][0] + arrow0[0]], [item.pt[1][1], item.pt[1][1] + arrow0[1]],
                             color='black')
                if np.linalg.norm(item.v[1]) >step_threshold and item.w[1] > step_threshold:
                    arrow1 = self.normalize(item.v[1]) *dt
                    plt.plot([item.pt[1][0], item.pt[1][0] + arrow1[0]], [item.pt[1][1], item.pt[1][1] + arrow1[1]],
                             color='blue')


        # plt.plot([self.nest_location[0], self.food_location[0]],
        #          [self.nest_location[1], self.food_location[1]], 'r*')
        plt.plot([default_nest_location[0], default_food_location[0]],
                 [default_nest_location[1], default_food_location[1]], 'r*')

        plt.plot([self.ants.ants[ant_tag].nt[1][0] for ant_tag in self.ants.ants if
                  self.ants.ants[ant_tag].mode[1] == 0],
                 [self.ants.ants[ant_tag].nt[1][1] for ant_tag in self.ants.ants if
                  self.ants.ants[ant_tag].mode[1] == 0], 'g*', markersize=2)
        plt.plot([self.ants.ants[ant_tag].nt[1][0] for ant_tag in self.ants.ants if
                  self.ants.ants[ant_tag].mode[1] == 1],
                 [self.ants.ants[ant_tag].nt[1][1] for ant_tag in self.ants.ants if
                  self.ants.ants[ant_tag].mode[1] == 1], 'y*', markersize=2)


        plt.xlim(-2, self.grid.domain[0]+2)
        plt.ylim(-2, self.grid.domain[1]+2)
        # plt.colorbar()
        if to_plot == 'W1' and fig_tag:
            plt.savefig(FOLDER_LOCATION + 'W1_WEIGHTS/' + str(to_plot) + '_' + str(fig_tag) + '.png')
            plt.close()
        elif to_plot == 'W2' and fig_tag:
            plt.savefig(FOLDER_LOCATION +'W2_WEIGHTS/' + str(to_plot) + '_' + str(fig_tag) + '.png')
            plt.close()
        elif to_plot == 'W' and fig_tag:
            plt.savefig(FOLDER_LOCATION +'W_WEIGHTS/' + str(to_plot) + '_' + str(fig_tag) + '.png')
            plt.close()
        else:
            plt.show()
            plt.close()


    def store_nr_trips(self,t):
        # self.total_trips[t] = sum([self.ants.ants[ant_tag].trips for ant_tag in self.ants.ants])
        if t >0:
            self.total_trips[t] = max(sum([self.ants.ants[ant_tag].trips for ant_tag in self.ants.ants]), self.total_trips[t-1])
        else:
            self.total_trips[t] = sum([self.ants.ants[ant_tag].trips for ant_tag in self.ants.ants])

    def plot_trips(self,total_time,fig_tag=None):
        trips_sequence = np.array([self.total_trips[time] for time in range(0,total_time)]) / self.N_total

        plt.plot(np.array(range(0,total_time))*dt, trips_sequence, 'r')
        plt.xlabel("Time")
        plt.ylabel("#Trips / #Agents")

        if fig_tag:
            plt.savefig(FOLDER_LOCATION + 'total_trips_' + str(fig_tag) + '.png')
            plt.close()
        else:
            plt.show()



    def plt_range_beacons(self, fig_tag=None):
        # vor = Voronoi([item.pt[1] for item in self.beacons.beacons])
        fig, ax = plt.subplots(figsize=(12, 6))

        # if len(self.beacons.beacons) >3:
        #     vor = Voronoi([self.beacons.beacons[beac_tag].pt[1] for beac_tag in self.beacons.beacons])
        #     voronoi_plot_2d(vor, show_vertices=False)

        for beac_tag in self.beacons.beacons:
            circle = plt.Circle(self.beacons.beacons[beac_tag].pt[1], clip_range , color='r',fill=False)
            ax.add_patch(circle)


        for count, ant_tag in enumerate(self.ants.ants):
            if count ==1:
                circle = plt.Circle(self.ants.ants[ant_tag].nt[1], clip_range , color='g',fill=False)
                ax.add_patch(circle)

        plt.plot([self.nest_location[0], self.food_location[0]],
                 [self.nest_location[1], self.food_location[1]], 'r*')

        plt.plot([self.ants.ants[ant_tag].nt[1][0] for ant_tag in self.ants.ants if
                  self.ants.ants[ant_tag].mode[1] == 0],
                 [self.ants.ants[ant_tag].nt[1][1] for ant_tag in self.ants.ants if
                  self.ants.ants[ant_tag].mode[1] == 0], 'g*')
        plt.plot([self.ants.ants[ant_tag].nt[1][0] for ant_tag in self.ants.ants if
                  self.ants.ants[ant_tag].mode[1] == 1],
                 [self.ants.ants[ant_tag].nt[1][1] for ant_tag in self.ants.ants if
                  self.ants.ants[ant_tag].mode[1] == 1], 'y*')

        plt.plot([self.beacons.beacons[beac_tag].pt[1][0] for beac_tag in self.beacons.beacons],
                 [self.beacons.beacons[beac_tag].pt[1][1] for beac_tag in self.beacons.beacons], 'b*')

        # plt.xlim(0, self.grid.domain[0])
        # plt.ylim(0, self.grid.domain[1])

        ax.set_xlim((0, self.grid.domain[0]))
        ax.set_ylim((0, self.grid.domain[1]))

        if fig_tag:
            plt.savefig(FOLDER_LOCATION + 'range_beacons' + str(fig_tag) + '.png')
            plt.close()
        else:
            plt.show()

    @staticmethod
    def normalize(item):
        return item / np.linalg.norm(item)