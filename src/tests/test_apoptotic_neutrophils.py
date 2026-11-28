import numpy as np

from grid import Grid
from environment import add_pro_inflammatory
from apoptotic_neutrophils import ApoptoticNeutrophil


def test_apoptotic_neutrophil_immobility():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create the inflammatory and anti-inflammatory fields
    c = np.zeros((100, 100))
    g = np.zeros((100, 100))

    # create an apoptotic neutrophil
    apoptotic_cell = ApoptoticNeutrophil(50, 50)

    # add the apoptotic neutrophil to the grid
    grid.add_agent(apoptotic_cell, 50, 50)

    # record its initial position
    initial_position = grid.get_location(apoptotic_cell)

    # update the apoptotic neutrophil
    apoptotic_cell.update(grid, c, 1)

    # retrieve its position after the update
    updated_position = grid.get_location(apoptotic_cell)

    # verify that it has not moved
    assert updated_position == initial_position

    print("apoptotic neutrophil immobility test passed")


def test_apoptotic_neutrophil_engulfment():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create an apoptotic neutrophil
    apoptotic_cell = ApoptoticNeutrophil(50, 50)

    # add the apoptotic neutrophil to the grid
    grid.add_agent(apoptotic_cell, 50, 50)

    # mark the apoptotic neutrophil as engulfed
    apoptotic_cell.get_engulfed()

    # remove the engulfed neutrophil from the grid
    grid.remove_agent(apoptotic_cell)

    # verify that the neutrophil is no longer on the grid
    assert apoptotic_cell not in grid.get_all_agents()

    # verify that its engulfed state was recorded
    assert apoptotic_cell.is_engulfed()

    print("apoptotic neutrophil engulfment test passed")


def test_apoptotic_neutrophil_necrosis():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create the inflammatory field
    c = np.zeros((100, 100))

    # create an apoptotic neutrophil
    apoptotic_cell = ApoptoticNeutrophil(50, 50)

    # add the apoptotic neutrophil to the grid
    grid.add_agent(apoptotic_cell, 50, 50)

    # set its lifespan to zero
    apoptotic_cell.lifespan = 0

    # trigger its update
    apoptotic_cell.update(grid, c, 1)

    # verify that the cell was removed from the grid
    assert apoptotic_cell not in grid.get_all_agents()

    # verify that its state changed to necrotic
    assert apoptotic_cell.state == "necrotic"

    print("apoptotic neutrophil necrosis test passed")


def test_apoptotic_neutrophil_mediator_release():

    # create a 100 x 100 grid
    grid = Grid(100, 100)

    # create the inflammatory field
    c = np.zeros((100, 100))

    # create an apoptotic neutrophil
    apoptotic_cell = ApoptoticNeutrophil(50, 50)

    # add the apoptotic neutrophil to the grid
    grid.add_agent(apoptotic_cell, 50, 50)

    # record the initial inflammatory concentration
    initial_concentration = c[50, 50]

    # trigger necrosis
    apoptotic_cell.undergo_necrosis(grid, c, 1)

    # retrieve the updated inflammatory concentration
    updated_concentration = c[50, 50]

    # verify that the inflammatory concentration increased
    assert updated_concentration > initial_concentration

    # verify that the increase equals the specified release amount
    assert np.isclose(updated_concentration, 1.0)

    print("apoptotic neutrophil mediator release test passed")