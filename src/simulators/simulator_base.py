import csv
import pathlib
import time
from abc import ABC
from enum import Enum

import matplotlib.path as mplt_path
import matplotlib.pyplot as plt
import numpy as np

from configs.base_simulation_config import BaseSimulationConfig
from utils.geometry import *
from utils.miscellaneous import *
from utils.scattering import *

from scipy.spatial import KDTree


class EdgeType(int, Enum):
    SourceEdge = 2,
    DrainEdge = 1
    MirrorEdge = 0
    RoughEdge = -1
    FalseEdge = -2


_SOURCE_EDGE = 2
_DRAIN_EDGE = 1
_MIRROR_EDGE = 0
_ROUGH_EDGE = -1
_FALSE_EDGE = -2


class SimulatorBase(ABC):
    """
    A base simulator class
    """

    def __init__(self, temperature: float, config: BaseSimulationConfig):
        self._temperature = temperature
        self._experimental_particle_density = config.experimental_particle_density
        self._simulation_particle_density = config.simulation_particle_density
        self._source_drain_ratio = config.source_drain_ratio
        self._e_min = calc_emin(self._experimental_particle_density,
                                self._simulation_particle_density,
                                self._temperature)

        self._p_scatter = config.scattering_probability

        self._overlap_radius = config.collision_distance

        # TODO: This value does not seem to be changed.
        # Could be smaller.
        self._initialization_count = 500_000

        # It is self.DeltaX in previous code
        self._step_size = config.step_size

        self._re_inject_with_unit_energy = config.re_inject_with_unit_energy

        # TODO: Does not seem to be changed
        self._number_drain_bins = 10
        self._n_corners_crossed = 0

        self._tip = config.add_probe_tip

        # TODO: Make this configurable
        self._probe_tip_diameter = 1.5
        self._probe_center_x = config.probe_center_x
        self._probe_center_y = config.probe_center_y

        self._rho_matrix = np.zeros([])
        self._px_matrix = np.zeros([])
        self._py_matrix = np.zeros([])
        self._e_rho_matrix = np.zeros([])

        self._n_restarts = 0
        self._time_count = 0

        self._generate_re_injection_probs = config.re_injection_probabilities
        self._probe_method = config.probe_method

        # TODO: Implement support for generated probability
        self._generated_probs = []
        self._roll_index = 0

        self.update_body()
        self.calc_area()
        self.set_field_resolution(config.field_resolution)
        self._set_n_part(self._simulation_particle_density, self._area)

        # self._overlaps = np.zeros((self._n_part, self._n_part), dtype=bool)

        self._update_particles()
        # TODO: This looks redundant, as it is already invoked in _update_particles()
        # self._find_overlaps()

        self._save_interval = config.save_interval

        self._make_movies = config.make_movie
        self._counts_per_snap_shot = 1000
        self._frame_num = 0

        # Keeps track of pairs that were overlapping (within collision distance) in previous time-step.
        self._overlapping_pairs = set()

    def _set_boundary_length(self, length):
        self._box_l = length
        self._box_range = cal_boundary_length(length)

    def _update_scattering_probability(self, value):
        self._p_scatter = value

    def calc_scatter_probability(self):
        # TODO: Where did 0.000024791 come from?
        self._p_scatter = 0.000024791 * self._temperature

    def update_body(self):
        # x-width and y-width???
        xw = self._box_l / 2.0
        yw = self._box_l / 2.0

        # Left border
        self._border_x = np.array([-xw, -xw])
        self._border_y = np.array([-yw, yw])

        # Create binned (???) drain edges (top)
        bin_length = self._box_l / self._number_drain_bins
        for bin_index in np.arange(1, self._number_drain_bins):
            self._border_x = np.append(self._border_x, self._border_x[bin_index] + bin_length)
            self._border_y = np.append(self._border_y, yw)

        self._border_x = np.append(self._border_x, [[xw, xw, -xw]])
        self._border_y = np.append(self._border_y, [[yw, -yw, -yw]])

        # One mirror edge + [a drain edge for each drain bin] + mirror edge + source_edge
        # self._edge_styles = [_MIRROR_EDGE] + [_DRAIN_EDGE] * self._number_drain_bins + [_MIRROR_EDGE, _SOURCE_EDGE]
        self._edge_styles = [_MIRROR_EDGE]
        for _ in np.arange(self._number_drain_bins):
            self._edge_styles += [_DRAIN_EDGE]
        self._edge_styles += [_MIRROR_EDGE, _SOURCE_EDGE]
        self._edge_styles = np.arange(self._edge_styles)

        self._device_edges_count = len(self._edge_styles)

        self._setup_edges()

        if self._tip:
            self._add_tip()

    def _setup_edges(self):
        corner_count = len(self._edge_styles)

        print(f"Corner count: {corner_count}")

        self._norm_x, self._norm_y, self._lengths = poly_norms(self._border_x, self._border_y)

        self._source_indices = [index for index, value in enumerate(self._edge_styles) if value == _SOURCE_EDGE]
        self._drain_indices = [index for index, value in enumerate(self._edge_styles) if value == _DRAIN_EDGE]

        self._probs = np.zeros(corner_count)

        if self._probe_method == "S":
            # TODO: Implement this after finding out about self.generatedProbs passed as argument to the program.
            raise NotImplementedError("Probe method by simulation is not supported.")
        else:
            source_length = np.sum(self._lengths[self._source_indices])
            drain_length = np.sum(self._lengths[self._drain_indices])

            # Calculate probability of injecting from source edges
            self._source_probs = (self._lengths[self._source_indices] * self._source_drain_ratio /
                                  (drain_length + source_length * self._source_drain_ratio))

            # Calculate probability of injecting from drain edges
            self._drain_probs = (self._lengths[self._drain_indices] /
                                 (drain_length + source_length * self._source_drain_ratio))

        self._probs[self._source_indices] = self._source_probs
        self._probs[self._drain_indices] = self._drain_probs

        # indices with > 0 probabilities
        self._prob_idx = (np.where(self._probs > 0))[0]

        # cumulative sum of probabilities, per index. i.e, prob[i] += prob[i-1]
        self._cumulative_prob = np.cumsum(self._probs[self._prob_idx])

        self._absorbed = np.zeros(corner_count)
        self._injected = np.zeros(corner_count)

        self._symmetrized_injected_total_old = np.zeros(self._device_edges_count)
        self._symmetrized_injected_new = np.zeros(self._device_edges_count)
        self._symmetrized_injected_old = np.zeros(self._device_edges_count)
        self._symmetrized_injected_dif = np.zeros(self._device_edges_count)

        # Number of ohmic contacts (M)
        contacts_count = len(self._source_indices) + len(self._drain_indices)

        # (M x M) matrix for 'Nabs' (number of absorbed??) between pairs of ohmic contacts;
        # rows = inj, columns = abs
        self._current_matrix = np.zeros((contacts_count, contacts_count))

        # sorted list of ohmic contact indices
        self._contact_lookup = np.sort(np.concatenate([self._source_indices, self._drain_indices]))

        self._border_path = mplt_path.Path(np.vstack((self._border_x, self._border_y)).T)

    def _add_tip(self):
        """
        Add circular probe tip
        :return:
        """

        # Equally spaced 50 numbers from pi -> 3pi
        thetas = np.linspace(np.pi, 3 * np.pi, 50)

        self._border_x = np.append(self._border_x, self._probe_center_x + np.cos(thetas) * self._probe_tip_diameter / 2)
        self._border_x = np.append(self._border_x, self._border_x[0])
        self._border_y = np.append(self._border_y, self._probe_center_y + np.sin(thetas) * self._probe_tip_diameter / 2)
        self._border_y = np.append(self._border_y, self._border_y[0])

        self._edge_styles = np.append(self._edge_styles, _FALSE_EDGE)
        for angle in thetas:
            if angle == 3 * np.pi:
                break
            self._edge_styles = np.append(self._edge_styles, _MIRROR_EDGE)
        self._edge_styles = np.append(self._edge_styles, _FALSE_EDGE)

    def calc_area(self):
        """
        Calculates the area of the device.
        :return:
        """
        false_edges_indices = np.where(self._edge_styles == _FALSE_EDGE)[0] + 1  # Adds 1 to each index
        border_x_split = np.split(self._border_x, false_edges_indices)  # Splits x-border at false edges
        border_y_split = np.split(self._border_y, false_edges_indices)  # Splits y-border at false edges

        area = 0
        for comp_index, x_component in enumerate(border_x_split):
            y_component = border_y_split[comp_index]
            for index in range(len(x_component)):
                next_index = index + 1 if index < len(x_component) - 1 else 0
                area += x_component[index] * y_component[next_index] - x_component[next_index] * y_component[index]

        self._area = np.abs(area / 2)
        print(f"Setting area: {self._area}")
        self._border_x_split = border_x_split
        self._border_y_split = border_y_split

    def set_field_resolution(self, dx: float):
        self._field_resolution = dx

        self._row_size = int(np.round(self._box_l / dx))
        self._column_size = int(np.round(self._box_l / dx))

        self._rho_matrix = np.zeros((self._row_size, self._column_size))
        self._px_matrix = np.zeros((self._row_size, self._column_size))
        self._py_matrix = np.zeros((self._row_size, self._column_size))
        self._e_rho_matrix = np.zeros((self._row_size, self._column_size))

        self._transport_map = np.array([np.zeros((self._row_size, self._column_size))] * len(self._contact_lookup))

        _, self._hist_X, self._hist_Y = np.histogram2d(
            np.zeros(self._row_size), np.zeros(self._column_size), bins=[self._row_size, self._column_size],
            range=self._box_range
        )

    def _set_n_part(self, simulating_particle_density, area):
        print(f"Density: {simulating_particle_density}, Area: {area}")
        self._n_part = int(simulating_particle_density * area)
        print(f"Total particles: {self._n_part}")
        self._injection_look_up = -1 * np.ones(self._n_part)

    def _update_particles(self):

        # Randomly assigns (x, y) position to N particles
        # TODO: Why -0.5?
        self._x_pos = self._box_l * (np.random.rand(self._n_part) - 0.5)
        self._y_pos = self._box_l * (np.random.rand(self._n_part) - 0.5)

        # out of bounds particles, points in (_x_pos, _y_pos) not within the area enclosed by _border_path
        print(f"Finding out of bound particles.")
        out_of_bounds = find_out_of_bounds(self._x_pos, self._y_pos, self._border_path)
        while len(out_of_bounds) > 0:
            for oob_index in out_of_bounds:
                self._x_pos[oob_index] = self._box_l * (np.random.rand() - 0.5)
                self._y_pos[oob_index] = self._box_l * (np.random.rand() - 0.5)
            out_of_bounds = find_out_of_bounds(self._x_pos, self._y_pos, self._border_path)

        self._overlapping_pairs = self._find_overlapping_pairs()

        # Assign random velocity direction
        thetas = np.random.rand(self._n_part) * 2.0 * np.pi
        self._v_x = np.cos(thetas)
        self._v_y = np.sin(thetas)

        # Assign fixed momentum amplitude
        self._p_r = np.ones(np.shape(thetas))

        self._l_no_scatter = np.zeros(self._n_part)
        self._traces = [[]] * self._n_part  # a list of #n_part_in empty lists

    def _find_overlapping_pairs(self):
        kd_tree = KDTree(np.vstack((self._x_pos, self._y_pos)).T)
        return kd_tree.query_pairs(self._overlap_radius)

    def _find_inject_position(self, contact_idx=None):
        """
        Calculates position to inject the given contact index. If contact index is missing, find one based
        on pre-calculated probability.
        :return:
        """
        # TODO: Is condition inside np.where() evaluated per item in the array?
        if not contact_idx:
            contact_idx = self._prob_idx[np.where(self._cumulative_prob > np.random.rand())[0][0]]
        inject_fraction = np.random.rand()

        x_1 = self._border_x[contact_idx + 1]
        x_0 = self._border_x[contact_idx]

        y_1 = self._border_y[contact_idx + 1]
        y_0 = self._border_y[contact_idx]

        inject_x = x_0 + inject_fraction * (x_1 - x_0)
        inject_y = y_0 + inject_fraction * (y_1 - y_0)

        return inject_x, inject_y, contact_idx

    def consume_and_re_inject(self, particle_index, current_edge_index):
        new_edge_index = current_edge_index if self._generate_re_injection_probs or \
                                               self._source_drain_ratio == 1.0 else None
        x_new, y_new, new_edge_index = self._find_inject_position(new_edge_index)

        self._absorbed[current_edge_index] += 1
        self._injected[new_edge_index] += 1

        self._update_current_matrix(particle_index, current_edge_index, new_edge_index)

        v_x_new, v_y_new = rand_cos_norm(self._norm_x[new_edge_index], self._norm_y[new_edge_index])

        self._v_x[particle_index] = v_x_new
        self._v_y[particle_index] = v_y_new

        # If re-injecting with unit energy.
        if self._re_inject_with_unit_energy:
            self._p_r[particle_index] = 1

        # Small offset to make sure the particle is in bounds.
        self._x_pos[particle_index] = x_new + v_x_new * self._step_size * 0.0001
        self._y_pos[particle_index] = y_new + v_y_new * self._step_size * 0.0001

    def _update_current_matrix(self, particle_index, absorbed_contact_index, injected_contact_index):
        """
        Tabulates (previous inject, new absorption) statistics
        :param particle_index:
        :param absorbed_contact_index:
        :param injected_contact_index:
        :return:
        """
        prev_injected_contact_index = self._injection_look_up[particle_index]
        if prev_injected_contact_index != -1:
            prev_injected_contact = np.where(self._contact_lookup == prev_injected_contact_index)
            prev_absorbed_contact = np.where(self._contact_lookup == absorbed_contact_index)
            self._current_matrix[prev_injected_contact, prev_absorbed_contact] += 1
        self._injection_look_up[particle_index] = injected_contact_index

    def _set_diffusive_edges(self):
        print(f"Setting diffusive edges")
        # Changes mirror edges to rough edges
        self._edge_styles = [_ROUGH_EDGE if edge == _MIRROR_EDGE else edge for edge in self._edge_styles]
        self._setup_edges()

    def _set_mirror_edges(self):
        print(f"Setting mirror edges")
        # Changes rough edges to mirror edges
        self._edge_styles = [_MIRROR_EDGE if edge == _ROUGH_EDGE else edge for edge in self._edge_styles]
        self._setup_edges()

    def _scatter_from_diffusive_edge(self, particle_id: int, crossed_edge_id: int):
        self._v_x[particle_id], self._v_y[particle_id] = rand_cos_norm(
            self._norm_x[crossed_edge_id], self._norm_y[crossed_edge_id]
        )

    def _reflect_from_mirror_edge(self, particle_id: int, crossed_edge_id: int):
        self._v_x[particle_id], self._v_y[particle_id] = mirror_norm(
            self._norm_x[crossed_edge_id], self._norm_y[crossed_edge_id],
            self._v_x[particle_id], self._v_y[particle_id]
        )

    def propagate_unit_time(self):
        # Propagate all particles in the direction of their velocity
        self._x_pos += self._v_x * self._step_size
        self._y_pos += self._v_y * self._step_size

    def backtrack(self, out_of_bound_indices):
        # Backtrack out of bound particles' position and capture new position
        oob_x = self._x_pos[out_of_bound_indices] - self._v_x[out_of_bound_indices] * self._step_size
        oob_y = self._y_pos[out_of_bound_indices] - self._v_y[out_of_bound_indices] * self._step_size

        # Recompute bound status of out of bound particles, using backtracked position of particles.
        # remaining_oob_indices = find_out_of_bounds(oob_x, oob_y, self._border_path)
        i, remaining_oob_indices = find_bound_status(oob_x, oob_y, self._border_x, self._border_y)
        while remaining_oob_indices.size > 0:
            # print(f"Remaining out of bounds: {len(remaining_oob_indices)}")
            thetas = np.random.rand(remaining_oob_indices.size) * 2 * np.pi
            oob_x[remaining_oob_indices] = \
                self._x_pos[out_of_bound_indices[remaining_oob_indices]] - np.sin(thetas) * self._step_size
            oob_y[remaining_oob_indices] = \
                self._y_pos[out_of_bound_indices[remaining_oob_indices]] - np.cos(thetas) * self._step_size

            # Ideally, we should only check for remaining points, instead of oob_x and oob_y. But then we will lose
            # their relative position in oob_x and oob_y, which is needed to backtrack them.
            # remaining_oob_indices = find_out_of_bounds(oob_x, oob_y, self._border_path)
            i, remaining_oob_indices = find_bound_status(oob_x, oob_y, self._border_x, self._border_y)
        return oob_x, oob_y

    @staticmethod
    def categorize_crossed_edges(crossed_edges_list, out_of_bound_indices):
        none_crossed = [out_of_bound_indices[val[0]] for val in crossed_edges_list if len(val[1]) < 1]
        one_crossed = [(out_of_bound_indices[val[0]], val[1][0]) for val in crossed_edges_list if len(val[1]) == 1]
        multiple_crossed = [(out_of_bound_indices[val[0]], val[1]) for val in crossed_edges_list if len(val[1]) > 1]
        return multiple_crossed, none_crossed, one_crossed

    def time_step(self):
        """
        Executes a single timestep of a simulation.
        """

        self.propagate_unit_time()

        # Update mean free path statistics
        self._l_no_scatter += self._step_size

        # Find particles (their indices) that are inside and outside the border path.
        # within_bound_indices, out_of_bound_indices = get_inbound_outbound(self._x_pos, self._y_pos, self._border_path)
        within_bound_indices, out_of_bound_indices = find_bound_status(self._x_pos, self._y_pos, self._border_x, self._border_y)

        # Step 1: Take care of in bounds particles.
        randomly_scatter(within_bound_indices, self._p_scatter, self._v_x, self._v_y)

        # Step 2: Take care of out of bound particles.

        oob_x, oob_y = self.backtrack(out_of_bound_indices)

        # It is defined here so that it can access the local variables (e.g., oob_x) defined above.
        # TODO: Make this operate on an array and Simplify
        def find_crossed_edges(index):
            # Returning index to avoid relying on order of result array.
            return [index, seg_cross_poly(self._border_x,
                                          self._border_y,
                                          np.array((oob_x[index], oob_y[index])),
                                          np.array((self._x_pos[out_of_bound_indices[index]],
                                                    self._y_pos[out_of_bound_indices[index]])))]

        # Make the function numpy uFunc (Universal Function) so that it can be executed element-wise.
        crossed_edges_list = self.perform_crossed_edges(find_crossed_edges, out_of_bound_indices)

        self._x_pos[out_of_bound_indices] = oob_x
        self._y_pos[out_of_bound_indices] = oob_y

        # Categorize the particles based on number of edges they crossed.
        multiple_crossed, none_crossed, one_crossed = \
            self.categorize_crossed_edges(crossed_edges_list, out_of_bound_indices)

        # Handle particles that did not cross any edge. Exactly between two.
        self._v_x[out_of_bound_indices[none_crossed]] *= -1
        self._v_y[out_of_bound_indices[none_crossed]] *= -1

        injected_particles = set()

        # Handle particles that crossed single edge.
        def scatter_particles_off_edge(ids):
            p_id = ids[0]
            edge_id = ids[1]
            # print(f"Edge id: {edge_id}")
            # print(f"Particle id: {p_id}")
            if self._edge_styles[edge_id] == _ROUGH_EDGE:  # Diffusive edge
                self._scatter_from_diffusive_edge(p_id, edge_id)
            if self._edge_styles[edge_id] == _MIRROR_EDGE:
                self._reflect_from_mirror_edge(p_id, edge_id)
            if self._edge_styles[edge_id] > 0:  # Source or Drain edges
                self.consume_and_re_inject(p_id, edge_id)
                injected_particles.add(p_id)

        if len(one_crossed) > 0:
            np.apply_along_axis(scatter_particles_off_edge, 1, one_crossed)

        self._n_corners_crossed += len(multiple_crossed)

        # Find an edge to scatter particles that crossed multiple edges.
        particle_edge_to_scatter_list = self.find_edge_for_multiple_crossed(multiple_crossed)

        if len(particle_edge_to_scatter_list) > 0:
            np.apply_along_axis(scatter_particles_off_edge, 1, particle_edge_to_scatter_list)

        overlapping_pairs = self._find_overlapping_pairs()

        for pair in overlapping_pairs:
            '''
            We do not process all overlapping pairs because the pairs can remain overlapped even after their collision
            as the distance travelled per time-step is less than overlapping bound (0.01 vs 0.05). So, we only scatter 
            pairs that overlap for the first time. self._overlapping_pairs keeps track of those overlapping pairs in 
            previous time-steps, hence the false check below to identify first time overlapping pairs.
            '''
            i = pair[0]
            j = pair[1]
            to_ignore = pair in self._overlapping_pairs or i in injected_particles or j in injected_particles

            if not to_ignore:
                self.scatter_particle_pair(i, j)

            self._overlapping_pairs.add((i, j))

        # Replace with the one just used for next time.
        self._overlapping_pairs = overlapping_pairs

        self._time_count += 1

        if self._time_count == self._initialization_count:
            self.reset_data()

        # Update histogram every 20th step
        if self._time_count % 20 == 0:
            self.collect_data()

    def reset_data(self):
        # Reset statistics after initialization count steps. TODO: Why?
        # TODO: May be cheaper to set new matrix instead of multiplying by 0
        self._px_matrix *= 0
        self._py_matrix *= 0
        self._rho_matrix *= 0
        self._e_rho_matrix *= 0
        self._injected *= 0
        self._absorbed *= 0
        self._current_matrix *= 0
        self._transport_map *= 0

    def collect_data(self):
        # TODO: box_range only includes [-L/2, L/2]. Why do we not consider points beyond this range?
        # TODO: What is seamless way to collect data, than calculating histogram?
        # TODO: Find what do we this for, and based on our usage, we might find cheaper ways to collect data.
        h, _, _ = np.histogram2d(self._x_pos, self._y_pos, bins=[self._row_size, self._column_size],
                                 range=self._box_range, weights=self._v_x)
        self._px_matrix += h
        h, _, _ = np.histogram2d(self._x_pos, self._y_pos, bins=[self._row_size, self._column_size],
                                 range=self._box_range, weights=self._v_y)
        self._py_matrix += h
        h, _, _ = np.histogram2d(self._x_pos, self._y_pos, bins=[self._row_size, self._column_size],
                                 range=self._box_range, weights=self._p_r)
        self._e_rho_matrix += h
        h, _, _ = np.histogram2d(self._x_pos, self._y_pos, bins=[self._row_size, self._column_size],
                                 range=self._box_range)
        self._rho_matrix += h
        # Update particle density for each Ohmic edge.
        # Only re-injected particles from that edge contribute.
        for index, edge in enumerate(self._contact_lookup):
            h, _, _ = np.histogram2d(self._x_pos[self._injection_look_up == edge],
                                     self._y_pos[self._injection_look_up == edge],
                                     bins=[self._row_size, self._column_size],
                                     range=self._box_range)
            self._transport_map[index] += h

    def scatter_particle_pair(self, i, j):
        # Perform relativistic scattering
        p1 = np.array([self._v_x[i] * self._p_r[i], self._v_y[i] * self._p_r[i]])
        p2 = np.array([self._v_x[j] * self._p_r[j], self._v_y[j] * self._p_r[j]])
        # Calculate outgoing momenta
        p3, p4 = scatter_massless_particles(p1, p2, self._e_min)
        if not (np.array_equal(p1, p3) and np.array_equal(p2, p4)):
            # Measure momentum amplitude
            self._p_r[i] = np.sqrt(np.sum(p3 ** 2))
            self._p_r[j] = np.sqrt(np.sum(p4 ** 2))

            # Set velocities
            self._v_x[i], self._v_y[i] = p3 / self._p_r[i]
            self._v_x[j], self._v_y[j] = p4 / self._p_r[j]

            self._l_no_scatter[i] = 0
            self._l_no_scatter[j] = 0

    def find_edge_for_multiple_crossed(self, multiple_crossed):
        particle_edge_to_scatter_list = []
        for item in multiple_crossed:
            particle_id = item[0]
            crossed_edges = item[1]
            crossed_source_drains = np.where(self._edge_styles[crossed_edges] > 0)[0]
            if len(crossed_source_drains):
                random_edge = np.random.randint(len(crossed_source_drains))
                particle_edge_to_scatter_list.append((particle_id, random_edge))
            else:
                # Choose the first edge whose norm has a component in the direction of velocity
                index = 0
                selected_edge = crossed_edges[index]
                while (self._v_x[particle_id] * self._norm_x[selected_edge] +
                       self._v_y[particle_id] * self._norm_y[selected_edge]) < 0:
                    index += 1
                    selected_edge = crossed_edges[index]
                particle_edge_to_scatter_list.append((particle_id, selected_edge))
        return particle_edge_to_scatter_list

    def perform_crossed_edges(self, find_crossed_edges, out_of_bound_indices):
        find_crossed_edges_func = np.frompyfunc(find_crossed_edges, 1, 1)
        crossed_edges_list = find_crossed_edges_func(np.arange(out_of_bound_indices.size))
        # Filter false edges from the crossed edges.
        for i in range(len(crossed_edges_list)):
            edges_to_keep = [crossed_edge for crossed_edge in crossed_edges_list[i][1]
                             if self._edge_styles[crossed_edge] != _FALSE_EDGE]
            crossed_edges_list[i][1] = edges_to_keep
        return crossed_edges_list

    # TODO: In old code filename is not being used. See how it is working.
    def _save_state(self, filename: pathlib.Path):
        filename.touch(exist_ok=True)
        var_names = list(self.__dict__.keys())
        save_func = f"np.savez('{filename}', " + ", ".join(map(lambda name: f"{name}=self.{name}", var_names)) + ")"
        exec(save_func)

    def _save_re_injection_stats(self, file_name: pathlib.Path):
        # This only applies to symmetric devices.
        # The roll index is unique to every geometry and can be specified in that device's init function
        injected_device_edges_only = np.roll(self._injected[0: self._device_edges_count], self._roll_index)

        self._symmetrized_injected_new = injected_device_edges_only - self._symmetrized_injected_total_old
        self._symmetrized_injected_dif = self._symmetrized_injected_new - self._symmetrized_injected_old

        with open(file_name, 'a') as file:
            csv_writer = csv.writer(file)
            csv_writer.writerows([self._symmetrized_injected_new.tolist()])

        self._symmetrized_injected_total_old = injected_device_edges_only
        self._symmetrized_injected_old = self._symmetrized_injected_new

    def _save_frame(self, output_dir: pathlib.Path):
        file_name = output_dir / f"frame_{self._frame_num}"
        total_absorbed = np.sum(self._absorbed)
        fig, ax = plt.subplots(1, 1, figsize=(12, 12))

        for i, x_component in enumerate(self._border_x_split):
            y_component = self._border_y_split[i]
            ax.plot(x_component, y_component, 'black')

        ax.plot(self._x_pos[0], self._y_pos[0], 'b.', markersize=25)
        ax.plot(self._x_pos[1:], self._y_pos[1:], 'r.', markersize=5)

        ax.text(self._box_l / 5, self._box_l, f"total absorbed: {total_absorbed}", fontsize=14)
        ax.text(self._box_l / 5, self._box_l + 1.25, f"timestep: {self._time_count}", fontsize=14)
        ax.text(-self._box_l / 2, self._box_l, f"temperature: {self._temperature}", fontsize=14)
        ax.text(-self._box_l / 2, self._box_l + 1.25, f"sim density: {self._simulation_particle_density}", fontsize=14)
        ax.axis('off')
        fig.savefig(file_name, dpi=250)
        plt.close()
        self._frame_num += 1

    def run_and_save(self, time_steps: int, output_dir: pathlib.Path):
        state_file = output_dir / "state.npz"
        re_injection_stat_file = output_dir / "re_injection_stat.csv"

        self._save_state(state_file)

        timer = time.time()
        self._save_frame(output_dir)
        for step in range(time_steps):
            self.time_step()

            if self._make_movies:
                if self._time_count % self._save_interval == 0:
                    print("Saving frame")
                    self._save_frame(output_dir)

            if self._time_count % self._save_interval == 0:
                print(f"Step: {self._time_count}. Time: {time.time() - timer}")
                timer = time.time()
                self._save_state(state_file)

                # Every 10th of initialization count.
                if self._time_count % (self._initialization_count / 10) == 0:
                    if self._generate_re_injection_probs:
                        # Store re-injection stat
                        self._save_re_injection_stats(re_injection_stat_file)
