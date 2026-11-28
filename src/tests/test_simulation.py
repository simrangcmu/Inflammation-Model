import numpy as np
from unittest.mock import patch

from simulation import Simulation
from active_neutrophils import Neutrophil
from apoptotic_neutrophils import ApoptoticNeutrophil
from macrophage import Macrophage


# define shared parameters for the tests

neutrophil_params = {
    "movesPerTick_n": 2,
    "n_k_grad": 80,
    "n_sigma": 1.0,
    "prob_neutroReleasePro": 0.5,
    "inc_neutroReleasePro": 0.001
}

macrophage_params = {
    "movesPerTick_m": 1,
    "prob_macroPhagoApop": 1.0,
    "prob_macroReleaseAnti": 0.8,
    "inc_macroReleaseAnti": 0.001,
    "m_k_grad": 80,
    "m_sigma": 1.0,
    "prob_macroLeave": 0.8,
    "thresh_macroLeave": 0.02
}


def create_test_simulation():
    """
    create an empty simulation for testing.
    """

    return Simulation(
        grid_size=100,
        neutrophil_positions=[],
        macrophage_positions=[],
        neutrophil_params=neutrophil_params,
        macrophage_params=macrophage_params
    )


def test_initialization():
    """
    check that the simulation starts with the correct initial conditions.
    """

    # create a new simulation
    simulation = create_test_simulation()

    # check that simulation time starts at zero
    assert simulation.time == 0

    # check that the grid has the correct dimensions
    assert simulation.grid_size == 100

    # check that no agents are initially present
    assert len(simulation.neutrophils) == 0
    assert len(simulation.macrophages) == 0
    assert len(simulation.grid.get_all_agents()) == 0

    # check that the chemical fields were initialized
    assert simulation.environment[0][50, 50] == 1.0
    assert np.all(simulation.environment[1] == 0.0)


def test_neutrophil_recruitment_thresholds():
    """
    check that neutrophils are recruited only at locations
    satisfying the pro-inflammatory and anti-inflammatory thresholds.
    """

    # create a new simulation
    simulation = create_test_simulation()

    # clear the initial damage field
    simulation.environment[0][:] = 0.0
    simulation.environment[1][:] = 0.0

    # create three test locations
    simulation.environment[0][10, 10] = 0.1
    simulation.environment[0][20, 20] = 0.04
    simulation.environment[0][30, 30] = 0.1

    # make the third location fail the anti-inflammatory threshold
    simulation.environment[1][30, 30] = 0.02

    # force recruitment to succeed at every eligible location
    with patch("simulation.random.random", return_value=0.0):

        simulation.recruit_neutrophils()

    # only the first location should have recruited a neutrophil
    recruited_neutrophils = simulation.grid.get_agents_by_type(Neutrophil)

    assert len(recruited_neutrophils) == 1
    assert simulation.grid.get_location(recruited_neutrophils[0]) == (10, 10)


def test_macrophage_recruitment_thresholds():
    """
    check that macrophages are recruited only where the
    pro-inflammatory concentration exceeds the recruitment threshold.
    """

    # create a new simulation
    simulation = create_test_simulation()

    # clear the initial pro-inflammatory field
    simulation.environment[0][:] = 0.0

    # create three test locations
    simulation.environment[0][10, 10] = 0.5
    simulation.environment[0][20, 20] = 0.4
    simulation.environment[0][30, 30] = 0.3

    # force recruitment to succeed at every eligible location
    with patch("simulation.random.random", return_value=0.0):

        simulation.recruit_macrophages()

    # only the first location should have recruited a macrophage
    recruited_macrophages = simulation.grid.get_agents_by_type(Macrophage)

    assert len(recruited_macrophages) == 1
    assert simulation.grid.get_location(recruited_macrophages[0]) == (10, 10)


def test_neutrophil_population_limit():
    """
    check that neutrophil recruitment stops at the maximum
    population of 4000 cells.
    """

    # create a new simulation
    simulation = create_test_simulation()

    # fill the entire field with eligible recruitment conditions
    simulation.environment[0][:] = 1.0
    simulation.environment[1][:] = 0.0

    # create 4000 neutrophils at one location
    for agent_number in range(4000):

        neutrophil = Neutrophil(
            simulation.grid,
            simulation.environment,
            simulation.neutrophil_params
        )

        simulation.grid.add_agent(neutrophil, 50, 50)
        simulation.neutrophils.append(neutrophil)

    # attempt to recruit additional neutrophils
    with patch("simulation.random.random", return_value=0.0):

        simulation.recruit_neutrophils()

    # check that the population has not exceeded 4000
    total_neutrophils = len(
        simulation.grid.get_agents_by_type(Neutrophil)
    ) + len(
        simulation.grid.get_agents_by_type(ApoptoticNeutrophil)
    )

    assert total_neutrophils == 4000


def test_macrophage_population_limit():
    """
    check that macrophage recruitment stops at the maximum
    population of 1000 cells.
    """

    # create a new simulation
    simulation = create_test_simulation()

    # fill the entire field with eligible recruitment conditions
    simulation.environment[0][:] = 1.0

    # create 1000 macrophages at one location
    for agent_number in range(1000):

        macrophage = Macrophage(
            simulation.grid,
            simulation.environment,
            simulation.macrophage_params
        )

        simulation.grid.add_agent(macrophage, 50, 50)
        simulation.macrophages.append(macrophage)

    # attempt to recruit additional macrophages
    with patch("simulation.random.random", return_value=0.0):

        simulation.recruit_macrophages()

    # check that the population has not exceeded 1000
    total_macrophages = len(
        simulation.grid.get_agents_by_type(Macrophage)
    )

    assert total_macrophages == 1000


def test_recruitment_schedule():
    """
    check that neutrophil and macrophage recruitment are called
    at the correct simulation intervals.
    """

    # create a new simulation
    simulation = create_test_simulation()

    # replace recruitment methods with mock methods
    with patch.object(simulation, "recruit_neutrophils") as recruit_n, \
         patch.object(simulation, "recruit_macrophages") as recruit_m, \
         patch.object(simulation, "update_neutrophils"), \
         patch.object(simulation, "update_macrophages"), \
         patch.object(simulation, "update_apoptotic_neutrophils"), \
         patch.object(simulation, "update_environment"):

        # execute the first six simulation ticks
        for step_number in range(6):
            simulation.step()

    # neutrophil recruitment should occur at ticks 0, 2, and 4
    assert recruit_n.call_count == 3

    # macrophage recruitment should occur at ticks 0 and 5
    assert recruit_m.call_count == 2


def test_simulation_time_advances():
    """
    check that every simulation step advances time by one tick.
    """

    # create a new simulation
    simulation = create_test_simulation()

    # replace the update methods to isolate the time advancement
    with patch.object(simulation, "recruit_neutrophils"), \
         patch.object(simulation, "recruit_macrophages"), \
         patch.object(simulation, "update_neutrophils"), \
         patch.object(simulation, "update_macrophages"), \
         patch.object(simulation, "update_apoptotic_neutrophils"), \
         patch.object(simulation, "update_environment"):

        # advance the simulation by five ticks
        for step_number in range(5):
            simulation.step()

    # check that simulation time increased by five
    assert simulation.time == 5


def test_environment_updates():
    """
    check that the chemical fields are updated during a simulation step.
    """

    # create a new simulation
    simulation = create_test_simulation()

    # save the initial concentration
    initial_c = simulation.environment[0].copy()

    # replace agent updates and recruitment to isolate the PDE update
    with patch.object(simulation, "recruit_neutrophils"), \
         patch.object(simulation, "recruit_macrophages"), \
         patch.object(simulation, "update_neutrophils"), \
         patch.object(simulation, "update_macrophages"), \
         patch.object(simulation, "update_apoptotic_neutrophils"):

        # advance one simulation step
        simulation.step()

    # check that the pro-inflammatory field changed due to diffusion
    assert not np.array_equal(simulation.environment[0], initial_c)

    # check that the anti-inflammatory field remains zero
    assert np.all(simulation.environment[1] == 0.0)


def test_environment_array_references_are_preserved():
    """
    check that updating the environment does not replace the
    original concentration arrays.
    """

    # create a new simulation
    simulation = create_test_simulation()

    # save references to the original arrays
    original_c = simulation.environment[0]
    original_g = simulation.environment[1]

    # update the environment
    simulation.update_environment()

    # check that the original arrays are still being used
    assert simulation.environment[0] is original_c
    assert simulation.environment[1] is original_g
