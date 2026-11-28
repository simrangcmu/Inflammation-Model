import numpy as np

from environment import (
    create_environment,
    calculate_laplacian,
    update_pro_inflammatory,
    update_anti_inflammatory,
    update_fields,
    add_pro_inflammatory,
    add_anti_inflammatory,
    get_pro_inflammatory,
    get_anti_inflammatory,
)


def test_initial_conditions():
    """
    check that the initial mediator fields match the paper's
    specified initial conditions.
    """

    # initialize the environment
    c, g = create_environment()

    # check that both fields have the correct dimensions
    assert c.shape == (100, 100)
    assert g.shape == (100, 100)

    # check that the center of the damage region has c = 1
    assert c[50, 50] == 1.0

    # check that a point on the damage boundary is included
    assert c[50, 60] == 1.0

    # check that a point outside the damage region has c = 0
    assert c[50, 61] == 0.0

    # check that the anti-inflammatory field starts at zero everywhere
    assert np.all(g == 0.0)


def test_laplacian_calculation():
    """
    check that the five-point laplacian correctly identifies
    concentration differences between neighboring locations.
    """

    # create a field with zero concentration everywhere
    field = np.zeros((100, 100))

    # place one unit of concentration at the center
    field[50, 50] = 1.0

    # calculate the laplacian
    laplacian = calculate_laplacian(field)

    # the center should have a laplacian of -4
    assert laplacian[50, 50] == -4.0

    # the four directly neighboring locations should have a laplacian of 1
    assert laplacian[49, 50] == 1.0
    assert laplacian[51, 50] == 1.0
    assert laplacian[50, 49] == 1.0
    assert laplacian[50, 51] == 1.0

    # diagonal neighbors should have a laplacian of zero
    assert laplacian[49, 49] == 0.0


def test_periodic_boundaries():
    """
    check that the laplacian uses periodic boundaries, meaning
    that locations on opposite edges of the grid are neighbors.
    """

    # create a field with zero concentration everywhere
    field = np.zeros((100, 100))

    # place concentration at the top-left corner
    field[0, 0] = 1.0

    # calculate the laplacian
    laplacian = calculate_laplacian(field)

    # the location on the opposite horizontal edge is a neighbor
    assert laplacian[0, 99] == 1.0

    # the location on the opposite vertical edge is a neighbor
    assert laplacian[99, 0] == 1.0

    # check that coordinates wrap around to the opposite edge
    assert get_pro_inflammatory(field, 100, 0) == 1.0
    assert get_pro_inflammatory(field, 0, 100) == 1.0

    # check that negative coordinates also wrap correctly
    assert get_pro_inflammatory(field, -100, 0) == 1.0
    assert get_pro_inflammatory(field, 0, -100) == 1.0


def test_diffusion_and_decay():
    """
    check that diffusion and decay update both mediator fields
    using the explicit euler method.
    """

    # create uniform concentration fields
    c = np.full((100, 100), 2.0)
    g = np.full((100, 100), 3.0)

    # define numerical parameters
    dt = 0.1
    diffusion_c = 1.5
    diffusion_g = 1.5
    decay_c = 0.002
    decay_g = 0.002

    # update both fields
    updated_c, updated_g = update_fields(
        c, g,
        dt,
        diffusion_c, diffusion_g,
        decay_c, decay_g
    )

    # because the fields are uniform, their laplacians should be zero
    # therefore, only decay should change their concentrations

    expected_c = 2.0 * (1 - dt * decay_c)
    expected_g = 3.0 * (1 - dt * decay_g)

    # check the updated concentrations
    assert np.allclose(updated_c, expected_c)
    assert np.allclose(updated_g, expected_g)

    # check that the original fields were not modified
    assert np.all(c == 2.0)
    assert np.all(g == 3.0)


def test_mediator_additions():
    """
    check that pro-inflammatory and anti-inflammatory mediators
    can be added to the correct grid locations.
    """

    # create empty concentration fields
    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    # add mediator concentrations at the center
    add_pro_inflammatory(c, 50, 50, 0.5)
    add_anti_inflammatory(g, 50, 50, 0.3)

    # check that the concentrations increased correctly
    assert get_pro_inflammatory(c, 50, 50) == 0.5
    assert get_anti_inflammatory(g, 50, 50) == 0.3

    # add more mediator concentration at the same location
    add_pro_inflammatory(c, 50, 50, 0.2)
    add_anti_inflammatory(g, 50, 50, 0.2)

    # check that the new amounts were added to the existing concentrations
    assert np.isclose(get_pro_inflammatory(c, 50, 50), 0.7)
    assert np.isclose(get_anti_inflammatory(g, 50, 50), 0.5)

    # check that adding concentration at one location does not affect another
    assert get_pro_inflammatory(c, 49, 50) == 0.0
    assert get_anti_inflammatory(g, 49, 50) == 0.0