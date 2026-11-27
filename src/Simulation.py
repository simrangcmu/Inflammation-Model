import random
import numpy as np

from grid import Grid
from environment import create_environment, update_fields
from active_neutrophils import Neutrophil
from apoptotic_neutrophils import ApoptoticNeutrophil
from macrophage import Macrophage


class Simulation:

    def __init__(
        self,
        grid_size,
        neutrophil_positions,
        macrophage_positions,
        neutrophil_params,
        macrophage_params
    ):
        # initialize simulation time
        self.time = 0

        # store the grid dimensions
        self.grid_size = grid_size

        # initialize the chemical environment
        self.environment = create_environment()

        # create the spatial grid
        self.grid = Grid(grid_size, grid_size)

        # store the agent parameters
        self.neutrophil_params = neutrophil_params
        self.macrophage_params = macrophage_params

        # create lists to keep track of the agents
        self.neutrophils = []
        self.macrophages = []

        # initialize all agents
        self.initialize_neutrophils(neutrophil_positions)
        self.initialize_macrophages(macrophage_positions)

        self.dt = 0.01

        self.diffusion_c = 0.1
        self.diffusion_g = 0.1

        self.decay_c = 0.01
        self.decay_g = 0.01

    def step(self):
        """
        Advance the simulation by one time step.
        """

        # update all neutrophils
        self.update_neutrophils()

        # update all macrophages
        self.update_macrophages()

        # update the chemical environment
        self.update_environment()

        # advance simulation time
        self.time += 1

    def run(self, number_of_steps):
        """
        Run the simulation for a specified number of time steps.
        """

        for step_number in range(number_of_steps):

            # advance the simulation by one time step
            self.step()

            # print progress every 10 steps
            if self.time % 10 == 0:
                print("Simulation time:", self.time)
                print("Active neutrophils:", len(self.grid.get_agents_by_type(Neutrophil)))
                print("Macrophages:", len(self.grid.get_agents_by_type(Macrophage)))
                print("Apoptotic neutrophils:", ApoptoticNeutrophil.apoptotic_counter)

        activated_macrophages = sum(
        1 for macrophage in self.grid.get_agents_by_type(Macrophage)
        if macrophage.activated
    )

        print("Activated macrophages:", activated_macrophages)
        print("Macrophages:", len(self.grid.get_agents_by_type(Macrophage)))
                
    def initialize_neutrophils(self, positions):
        """
        Create neutrophils and place them at their initial positions.
        """

        # create each neutrophil at its specified position
        for x, y in positions:

            # create a new neutrophil
            neutrophil = Neutrophil(
                self.grid,
                self.environment,
                self.neutrophil_params
            )

            # place the neutrophil on the grid
            self.grid.add_agent(neutrophil, x, y)

            # keep track of the neutrophil
            self.neutrophils.append(neutrophil)

    def initialize_macrophages(self, positions):
        """
        Create macrophages and place them at their initial positions.
        """

        # create each macrophage at its specified position
        for x, y in positions:

            # create a new macrophage
            macrophage = Macrophage(
                self.grid,
                self.environment,
                self.macrophage_params
            )

            # place the macrophage on the grid
            self.grid.add_agent(macrophage, x, y)

            # keep track of the macrophage
            self.macrophages.append(macrophage)


    def update_neutrophils(self):
        """
        Update the behavior of all active neutrophils.
        """

        # create a copy of the current neutrophil list
        current_neutrophils = self.neutrophils.copy()

        # update each neutrophil
        for neutrophil in current_neutrophils:

            # only update neutrophils that are still on the grid
            if neutrophil in self.grid.get_all_agents():

                # run the neutrophil's behavior
                neutrophil.run()

        # remove neutrophils that are no longer on the grid
        self.neutrophils = [
            neutrophil
            for neutrophil in self.neutrophils
            if neutrophil in self.grid.get_all_agents()
    ]

    def update_macrophages(self):
        """
        Update the behavior of all macrophages.
        """

        # create a copy of the macrophage list
        current_macrophages = self.macrophages.copy()

        # update each macrophage
        for macrophage in current_macrophages:

            # only update macrophages that are still alive
            if macrophage.alive:
                macrophage.update()

    def update_environment(self):
        """
        Update the chemical fields using diffusion and decay.
        """

        updated_c, updated_g = update_fields(
            self.environment[0],
            self.environment[1],
            self.dt,
            self.diffusion_c,
            self.diffusion_g,
            self.decay_c,
            self.decay_g
        )

        # update the existing arrays so all agents retain the correct references
        self.environment[0][:] = updated_c
        self.environment[1][:] = updated_g

if __name__ == "__main__":

    # define temporary neutrophil parameters
    neutrophil_params = {
        "movesPerTick_n": 2,
        "n_k_grad": 0.5,
        "n_sigma": 1.0,
        "prob_neutroReleasePro": 0.5,
        "inc_neutroReleasePro": 0.1
    }

    # define temporary macrophage parameters
    macrophage_params = {
        "movesPerTick_m": 2,
        "prob_macroPhagoApop": 1.0,
        "threshold_macroActivation": 2,
        "prob_macroReleaseAnti": 1.0,
        "inc_macroReleaseAnti": 0.5,
        "m_k_grad": 0.5,
        "thresh_macroLeave": 0.5,
        "prob_macroLeave": 0.0
    }

    # define initial agent positions
    neutrophil_positions = [
        (45, 45),
        (48, 50),
        (55, 52),
        (52, 46)
    ]

    macrophage_positions = [
        (40, 40),
        (60, 60)
    ]

    # create the simulation
    simulation = Simulation(
        grid_size=100,
        neutrophil_positions=neutrophil_positions,
        macrophage_positions=macrophage_positions,
        neutrophil_params=neutrophil_params,
        macrophage_params=macrophage_params
    )

    # display initial conditions
    print("SIMULATION INITIALIZED")
    print("Active neutrophils:", len(simulation.neutrophils))
    print("Macrophages:", len(simulation.macrophages))

    # run the simulation for 100 time steps
    print("\nRUNNING SIMULATION")

    simulation.run(1500)

    # display final results
    print("\nSIMULATION COMPLETE")
    print("Final simulation time:", simulation.time)
    print("Active neutrophils:", len(simulation.neutrophils))
    print("Macrophages:", len(simulation.grid.get_agents_by_type(Macrophage)))
    print("Apoptotic neutrophils:", ApoptoticNeutrophil.apoptotic_counter)



    print("\nMACROPHAGE DIAGNOSTICS")

    for macrophage in simulation.grid.get_agents_by_type(Macrophage):

        print("\nMacrophage position:", simulation.grid.get_location(macrophage))
        print("Phagocytosis count:", macrophage.phagocytosis_count)
        print("Activated:", macrophage.activated)
        print("Alive:", macrophage.alive)


    print("\nAPOPTOTIC NEUTROPHIL POSITIONS")

    for apoptotic_neutrophil in simulation.grid.get_agents_by_type(ApoptoticNeutrophil):

        print(simulation.grid.get_location(apoptotic_neutrophil))


    print("\nCHEMICAL CONCENTRATIONS")

    print("Maximum pro-inflammatory concentration:", simulation.environment[0].max())
    print("Maximum anti-inflammatory concentration:", simulation.environment[1].max())