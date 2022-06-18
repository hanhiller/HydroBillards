import numpy as np
from numba import jit, types as t


@jit(t.Tuple((t.float64[:], t.float64))(t.float64[:], t.float64, t.float64[:]), nopython=True)
def _boost(p_bar, energy, v_bar):
    """
    Performs Lorentz Boost of momentum (p_bar) and energy (energy) along velocity (v_bar). (Velocity are considered in
    unit of 'c'.)

    Returns a new momentum and energy.
    :param p_bar:
    :param energy:
    :param v_bar:
    :return:
    """

    v_square = np.sum(v_bar ** 2)
    gamma = (1.0 - v_square) ** (-0.5)
    p_dot_v = np.sum(p_bar * v_bar)

    e_out = energy
    p_out = p_bar

    # Perform boost
    if v_square > 0:
        e_out = gamma * (energy - p_dot_v)
        p_out = p_bar + (gamma - 1.0) * p_dot_v * v_bar / v_square - gamma * energy * v_bar

    return p_out, e_out


@jit(t.float64[:](t.float64[:], t.float64[:]), nopython=True)
def _find_collinear_boost(p_bar_1, p_bar_2):
    """
    Computes a velocity vector needed for Lorentz Boost to make the two incoming momentum vectors collinear. It assumes
    particles are massless: and velocity is 'c' along the momentum direction.
    :param p_bar_1:
    :param p_bar_2:
    :return:
    """
    # Calculate direction of momentum
    theta1 = np.arctan2(p_bar_1[1], p_bar_1[0])
    theta2 = np.arctan2(p_bar_2[1], p_bar_2[0])

    # Find direction along which both vectors have the same projected amplitude.
    v_theta = np.arctan2((np.cos(theta1) - np.cos(theta2)), (np.sin(theta2) - np.sin(theta1)))

    # Calculate amplitude of the projected amplitude (???)
    v_amp = np.cos(theta1) * np.cos(v_theta) + np.sin(theta1) * np.sin(v_theta)

    # Output velocity vector to boost along
    return np.array([v_amp * np.cos(v_theta), v_amp * np.sin(v_theta)])


@jit(t.Tuple((t.float64[:], t.float64, t.float64[:], t.float64, t.float64[:], t.float64[:]))(t.float64[:], t.float64,
                                                                                             t.float64[:], t.float64),
     nopython=True)
def _find_center_of_mass_boost(p_bar_1, e1, p_bar_2, e2):
    """
    Takes momenta and energy of two particles and boosts them into a center of mass (COM) frame. It returns momenta
    and energy as well as two sequential boost velocities to get to the COM frame. The first velocity defines the
    boost required to get the vectors to be collinear and the second velocity is the boost required to get both
    particles to have the same energy.

    It currently assumes massless particles, since the COM frame is generally defined in the case of equal and
    opposite momenta.

    :param p_bar_1:
    :param e1:
    :param p_bar_2:
    :param e2:
    :return:
    """
    # Boost to bring momenta to be collinear
    v_boost_1 = _find_collinear_boost(p_bar_1, p_bar_2)
    p_bar_1_mid, e1_mid = _boost(p_bar_1, e1, v_boost_1)
    p_bar_2_mid, e2_mid = _boost(p_bar_2, e2, v_boost_1)

    # Calculate boost to get equal energies (and thus momenta for massless particles)
    p1 = np.sqrt(np.sum(p_bar_1_mid ** 2))
    p2 = np.sqrt(np.sum(p_bar_2_mid ** 2))
    v_boost_2 = p_bar_1_mid / p1 * (e1_mid - e2_mid) / (p1 + p2)

    # Perform the second boost
    p_bar_1_out, e1_out = _boost(p_bar_1_mid, e1_mid, v_boost_2)
    p_bar_2_out, e2_out = _boost(p_bar_2_mid, e2_mid, v_boost_2)

    return p_bar_1_out, e1_out, p_bar_2_out, e2_out, v_boost_1, v_boost_2


@jit(t.UniTuple(t.float64[:], 2)(t.float64[:], t.float64[:], t.float64), nopython=True)
# Relativistic Scattering
def scatter_massless_particles(p1, p2, e_min):
    e1 = np.sqrt(np.sum(p1 ** 2))
    e2 = np.sqrt(np.sum(p2 ** 2))

    # Boost to get to COM frame
    p10, e10, p20, e20, v1, v2 = _find_center_of_mass_boost(p1, e1, p2, e2)

    # Randomize output angles
    theta_0 = np.random.rand() * 2.0 * np.pi
    c0 = np.cos(theta_0)
    s0 = np.sin(theta_0)

    # Rotation matrix for rotating axes??
    rotation_matrix = np.array([[c0, s0], [-s0, c0]])

    p10 = p10 @ rotation_matrix
    # p10 = np.matmul(p10, rotation_matrix)
    p20 = p20 @ rotation_matrix
    # p20 = np.matmul(p20, rotation_matrix)

    # Reverse the second boost
    p11, e11 = _boost(p10, e10, -v2)
    p21, e21 = _boost(p20, e20, -v2)

    # Reverse the first boost
    p12, e12 = _boost(p11, e11, -v1)
    p22, e22 = _boost(p21, e21, -v1)

    # No scattering if energy less than e_min
    if e12 < e_min or e22 < e_min:
        return p1, p2
    else:
        if np.random.rand() < 0.5:
            return p12, p22
        else:
            return p22, p12


@jit(t.int64[:](t.int64[:], t.float64, t.float64[:], t.float64[:]), nopython=True)
def randomly_scatter(particle_ids, scattering_probability, v_x, v_y):
    # Randomly select particles to scatter, using scattering probability.
    indices_to_scatter = particle_ids[np.random.rand(particle_ids.size) <= scattering_probability]

    # Scatter selected particles
    thetas = np.random.rand(indices_to_scatter.size) * 2.0 * np.pi
    v_x[indices_to_scatter] = np.cos(thetas)
    v_y[indices_to_scatter] = np.sin(thetas)

    return indices_to_scatter
