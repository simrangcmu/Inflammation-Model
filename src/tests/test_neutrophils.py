import numpy as np

from grid import Grid
from active_neutrophils import Neutrophil, NEUTROPHIL_PARAMS
from apoptotic_neutrophils import ApoptoticNeutrophil


def test_neutrophil_chemotaxis():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create the inflammatory and anti-inflammatory fields
    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    # combine the two fields into the environment
    environment = (c, g)

    # create a neutrophil and place it at the center
    neutrophil = Neutrophil(grid, environment, NEUTROPHIL_PARAMS)
    grid.add_agent(neutrophil, 50, 50)

    # define the current position of the neutrophil
    x, y = grid.get_location(neutrophil)

    # assign an inflammatory concentration at the current location
    c[y, x] = 0.1

    # retrieve the eight neighboring locations
    neighbors = grid.get_neighbors(x, y)

    # assign different inflammatory concentrations to the neighbors
    # the right neighbor has the highest concentration
    c[50, 51] = 0.5

    # the upper-right neighbor has a smaller concentration
    c[51, 51] = 0.2

    # all remaining neighbors have a concentration of zero

    # calculate the probability of moving toward each neighbor
    neighbors, gradient_probabilities = neutrophil.calculate_gradient_probabilities()

    # check that all eight gradient probabilities add up to 1
    assert np.isclose(np.sum(gradient_probabilities), 1.0)

    # check that there are eight movement probabilities
    assert len(gradient_probabilities) == 8

    # identify the neighboring direction with the highest probability
    most_likely_direction = np.argmax(gradient_probabilities)

    # the right neighbor is the first direction in the directions list
    assert most_likely_direction == 0

    # check that moving right is more likely than moving upper-right
    assert gradient_probabilities[0] > gradient_probabilities[1]

    # check that the right neighbor has a higher probability than all others
    assert gradient_probabilities[0] > np.max(gradient_probabilities[1:])

    # calculate the combined movement probabilities
    movement_neighbors, movement_probabilities = neutrophil.calculate_movement_probabilities()

    # check that the movement probabilities sum to one
    assert np.isclose(np.sum(movement_probabilities), 1.0)

    # check that movement toward the right remains most likely
    assert np.argmax(movement_probabilities) == 0

    # verify that the initial memory probabilities are uniform
    memory_probabilities = neutrophil.calculate_memory_probabilities()

    memory_neighbors, memory_probabilities = neutrophil.calculate_memory_probabilities()

    # verify that the initial memory probabilities are uniform
    assert np.allclose(
        memory_probabilities,
        np.ones(8) / 8
    )

    print("neutrophil chemotaxis test passed")




def test_neutrophil_directional_persistence():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create uniform inflammatory and anti-inflammatory fields
    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    # combine the two fields into the environment
    environment = (c, g)

    # create a neutrophil and place it at the center
    neutrophil = Neutrophil(grid, environment, NEUTROPHIL_PARAMS)
    grid.add_agent(neutrophil, 50, 50)

    # set the previous movement direction to the right
    neutrophil.prev_move = 0

    # assign equal inflammatory concentrations to all eight neighbors
    neighbors = grid.get_neighbors(50, 50)

    for x, y in neighbors:
        c[y, x] = 0.1

    # calculate the memory probabilities
    memory_neighbors, memory_probabilities = neutrophil.calculate_memory_probabilities()

    # check that all eight probabilities add up to one
    assert np.isclose(np.sum(memory_probabilities), 1.0)

    # check that moving right has the highest memory probability
    assert np.argmax(memory_probabilities) == 0

    # check that continuing right is more likely than moving left
    assert memory_probabilities[0] > memory_probabilities[4]

    # calculate the combined movement probabilities
    movement_neighbors, movement_probabilities = neutrophil.calculate_movement_probabilities()

    # check that all movement probabilities add up to one
    assert np.isclose(np.sum(movement_probabilities), 1.0)

    # check that directional persistence favors moving right
    assert np.argmax(movement_probabilities) == 0

    print("neutrophil directional persistence test passed")


def test_neutrophil_moves_toward_highest_concentration():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create the inflammatory and anti-inflammatory fields
    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    # combine the fields into the environment
    environment = (c, g)

    # create a neutrophil at the center
    neutrophil = Neutrophil(grid, environment, NEUTROPHIL_PARAMS)
    grid.add_agent(neutrophil, 50, 50)

    # assign a high inflammatory concentration to the right
    c[50, 51] = 1.0

    # repeatedly test movement toward the gradient
    for _ in range(20):

        # reset the neutrophil to the center
        grid.move_agent(neutrophil, 50, 50)

        # move the neutrophil chemotactically
        neutrophil.move_chemotactically()

        # retrieve its new position
        new_x, new_y = grid.get_location(neutrophil)

        # verify that it moved to the right
        assert (new_x, new_y) == (51, 50)

    print("neutrophil movement test passed")


def test_neutrophil_becomes_apoptotic():

    # create a grid and environment
    grid = Grid(100, 100)

    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    environment = (c, g)

    # create an active neutrophil at the center
    neutrophil = Neutrophil(grid, environment, NEUTROPHIL_PARAMS)
    grid.add_agent(neutrophil, 50, 50)

    # assign a lifespan of zero to trigger the transition
    neutrophil.lifespan = 0

    # run the neutrophil update
    neutrophil.run()

    # retrieve the agents at the original location
    agents_at_location = grid.get_agents_at(50, 50)

    # verify that the active neutrophil was removed
    assert neutrophil not in grid.get_all_agents()

    # verify that an apoptotic neutrophil now occupies the location
    apoptotic_cells = [
        agent for agent in agents_at_location
        if isinstance(agent, ApoptoticNeutrophil)
    ]

    assert len(apoptotic_cells) == 1

    # verify that the apoptotic neutrophil remains at the same location
    assert grid.get_location(apoptotic_cells[0]) == (50, 50)

    print("neutrophil apoptosis transition test passed")

def test_neutrophil_pro_inflammatory_release():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create the inflammatory and anti-inflammatory fields
    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    # combine the fields into the environment
    environment = (c, g)

    # create a copy of the neutrophil parameters
    parameters = NEUTROPHIL_PARAMS.copy()

    # make mediator release certain for this test
    parameters["prob_neutroReleasePro"] = 1.0

    # create a neutrophil and place it at the center
    neutrophil = Neutrophil(grid, environment, parameters)
    grid.add_agent(neutrophil, 50, 50)

    # record the initial inflammatory concentration
    initial_concentration = c[50, 50]

    # trigger pro-inflammatory mediator release
    neutrophil.release_pro()

    # retrieve the updated concentration
    updated_concentration = c[50, 50]

    # verify that the concentration increased
    assert updated_concentration > initial_concentration

    print("neutrophil mediator release test passed")

def test_neutrophil_lifespan_decreases():

    # create a grid and environment
    grid = Grid(100, 100)

    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    environment = (c, g)

    # create a neutrophil
    neutrophil = Neutrophil(grid, environment, NEUTROPHIL_PARAMS)
    grid.add_agent(neutrophil, 50, 50)

    # assign a known lifespan
    neutrophil.lifespan = 10

    # decrease the lifespan
    neutrophil.update_lifespan()

    # verify that the lifespan decreased by one
    assert neutrophil.lifespan == 9

    print("neutrophil lifespan test passed")