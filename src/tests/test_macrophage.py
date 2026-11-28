import numpy as np

from grid import Grid
from macrophage import Macrophage, macrophage_params
from apoptotic_neutrophils import ApoptoticNeutrophil


def test_macrophage_chemotaxis():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create the inflammatory and anti-inflammatory fields
    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    # combine the fields into the environment
    environment = (c, g)

    # create a macrophage and place it at the center
    macrophage = Macrophage(grid, environment, macrophage_params.copy())
    grid.add_agent(macrophage, 50, 50)

    # create a strong inflammatory gradient toward the right
    c[50, 51] = 1.0

    # move the macrophage chemotactically
    macrophage.move_chemotactically()

    # retrieve the macrophage's new position
    new_position = grid.get_location(macrophage)

    # verify that it moved toward the right
    assert new_position == (51, 50)

    print("macrophage chemotaxis test passed")


def test_macrophage_moves_toward_apoptotic_neutrophil():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create the inflammatory and anti-inflammatory fields
    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    environment = (c, g)

    # create a macrophage at the center
    macrophage = Macrophage(grid, environment, macrophage_params.copy())
    grid.add_agent(macrophage, 50, 50)

    # create an apoptotic neutrophil to the right
    apoptotic_cell = ApoptoticNeutrophil(51, 50)
    grid.add_agent(apoptotic_cell, 51, 50)

    # move the macrophage
    macrophage.move()

    # retrieve its new position
    new_position = grid.get_location(macrophage)

    # verify that it moved toward the apoptotic neutrophil
    assert new_position == (51, 50)

    print("macrophage movement toward apoptotic neutrophil test passed")

def test_macrophage_phagocytosis():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create the inflammatory and anti-inflammatory fields
    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    environment = (c, g)

    # create a macrophage
    macrophage = Macrophage(grid, environment, macrophage_params.copy())
    grid.add_agent(macrophage, 50, 50)

    # create an apoptotic neutrophil at the same location
    apoptotic_cell = ApoptoticNeutrophil(50, 50)
    grid.add_agent(apoptotic_cell, 50, 50)

    # attempt phagocytosis
    macrophage.phagocytose()

    # verify that the apoptotic neutrophil was removed
    assert apoptotic_cell not in grid.get_all_agents()

    # verify that the macrophage successfully phagocytosed a cell
    assert macrophage.phagocytosis_count == 1

    print("macrophage phagocytosis test passed")

def test_macrophage_activation():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create the inflammatory and anti-inflammatory fields
    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    environment = (c, g)

    # create a macrophage
    macrophage = Macrophage(grid, environment, macrophage_params.copy())
    grid.add_agent(macrophage, 50, 50)

    # verify that the macrophage is initially inactive
    assert macrophage.activated is False

    # create an apoptotic neutrophil at the macrophage's location
    apoptotic_cell = ApoptoticNeutrophil(50, 50)
    grid.add_agent(apoptotic_cell, 50, 50)

    # perform successful phagocytosis
    macrophage.phagocytose()

    # verify that the macrophage is now activated
    assert macrophage.activated is True

    print("macrophage activation test passed")

def test_macrophage_anti_inflammatory_release():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create the inflammatory and anti-inflammatory fields
    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    environment = (c, g)

    # create a copy of the macrophage parameters
    parameters = macrophage_params.copy()

    # make anti-inflammatory mediator release certain
    parameters["prob_macroReleaseAnti"] = 1.0

    # create an activated macrophage
    macrophage = Macrophage(grid, environment, parameters)
    grid.add_agent(macrophage, 50, 50)

    # activate the macrophage
    macrophage.activated = True

    # record the initial anti-inflammatory concentration
    initial_concentration = g[50, 50]

    # trigger anti-inflammatory mediator release
    macrophage.release_mediator()

    # retrieve the updated concentration
    updated_concentration = g[50, 50]

    # verify that the concentration increased
    assert updated_concentration > initial_concentration

    print("macrophage anti-inflammatory release test passed")

def test_macrophage_departure():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create fields with no inflammatory signal
    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    environment = (c, g)

    # create a copy of the macrophage parameters
    parameters = macrophage_params.copy()

    # make macrophage departure certain when eligible
    parameters["prob_macroLeave"] = 1.0

    # create a macrophage
    macrophage = Macrophage(grid, environment, parameters)
    grid.add_agent(macrophage, 50, 50)

    # check whether the macrophage should leave
    macrophage.check_departure()

    # verify that the macrophage was removed from the grid
    assert macrophage not in grid.get_all_agents()

    # verify that the macrophage is no longer alive
    assert macrophage.alive is False

    print("macrophage departure test passed")


