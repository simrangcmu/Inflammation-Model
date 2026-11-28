import numpy as np
import random
from environment import add_pro_inflammatory
from apoptotic_neutrophils import ApoptoticNeutrophil

# the eight possible movement directions, represented as (dx, dy)
DIRECTIONS = [
    (1, 0),    # right
    (1, 1),    # up/right
    (0, 1),    # up
    (-1, 1),   # up/left
    (-1, 0),   # left
    (-1, -1),  # down/left
    (0, -1),   # down
    (1, -1)    # down/right
]


# neutrophil-specific parameters
NEUTROPHIL_PARAMS = {
    "movesPerTick_n": 2,
    "n_k_grad": 80,
    "n_sigma": 1.0,
    "prob_neutroReleasePro": 0.5,
    "inc_neutroReleasePro": 0.001,
}

class Neutrophil:

    neutrophil_counter = 0

    def __init__(self, grid, environment, params):

        # store references to the grid and environment
        self.grid = grid
        self.environment = environment
        self.params = params

        # initialize neutrophil properties
        self.lifespan = np.random.randint(60, 1441)
        self.prev_move = None

        # increment the total neutrophil count
        Neutrophil.neutrophil_counter += 1

    def get_neighbors(self):
        """
        retrieve the eight neighboring positions around the neutrophil.
        """

        # get the neutrophil's current location
        x, y = self.grid.get_location(self)

        # retrieve the neighboring coordinates from the grid
        neighbors = self.grid.get_neighbors(x, y)

        return neighbors

    def calculate_gradient_probabilities(self):
        """
        calculate movement probabilities based on differences
        in pro-inflammatory mediator concentration.
        """

        # retrieve the neutrophil's current location
        x, y = self.grid.get_location(self)

        # retrieve the neighboring positions
        neighbors = self.get_neighbors()

        # retrieve the pro-inflammatory concentration at the current location
        current_c = self.environment[0][y, x]

        # initialize a list to store the weights
        gradient_weights = []

        # calculate the weight for each neighboring position
        for neighbor in neighbors:

            # retrieve the neighboring coordinates
            neighbor_x, neighbor_y = neighbor

            # retrieve the pro-inflammatory concentration at the neighboring location
            neighbor_c = self.environment[0][neighbor_y, neighbor_x]

            # calculate the difference in concentration
            difference = neighbor_c - current_c

            # calculate the concentration differences for all neighbors
            concentration_differences = []

            for neighbor in neighbors:

                # retrieve the neighboring coordinates
                neighbor_x, neighbor_y = neighbor

                # retrieve the neighboring concentration
                neighbor_c = self.environment[0][neighbor_y, neighbor_x]

                # calculate the concentration difference
                difference = neighbor_c - current_c

                concentration_differences.append(difference)

            # convert the differences into a numpy array
            concentration_differences = np.array(concentration_differences)

            # subtract the largest difference to prevent numerical overflow
            adjusted_differences = (
                concentration_differences - np.max(concentration_differences)
            )

            # calculate the gradient weights
            gradient_weights = np.exp(
                self.params["n_k_grad"] * adjusted_differences
            )

            # normalize the weights to obtain probabilities
            gradient_probabilities = (
                gradient_weights / np.sum(gradient_weights)
            )

        return neighbors, gradient_probabilities
    
    def calculate_memory_probabilities(self):
        """
        calculate movement probabilities based on directional memory.

        directions closer to the previous movement direction
        receive higher probabilities.
        """

        # retrieve the possible neighboring positions
        neighbors = self.get_neighbors()

        # initialize a list to store memory weights
        memory_weights = []

        # retrieve the current location
        x, y = self.grid.get_location(self)
        
        # if the neutrophil has not moved yet, give every direction equal memory weight
        if self.prev_move is None:
            memory_probabilities = np.ones(len(neighbors)) / len(neighbors)
            return neighbors, memory_probabilities
        
    
        # retrieve the previous movement direction
        previous_dx, previous_dy = DIRECTIONS[self.prev_move]

        # calculate the angle of the previous movement
        previous_angle = np.arctan2(previous_dy, previous_dx)


        # calculate the memory weight for each possible movement
        for neighbor in neighbors:

            # retrieve neighboring coordinates
            neighbor_x, neighbor_y = neighbor

            # calculate the movement direction
            dx = neighbor_x - x
            dy = neighbor_y - y

            # correct for periodic boundaries
            if dx > 1:
                dx = -1
            elif dx < -1:
                dx = 1

            if dy > 1:
                dy = -1
            elif dy < -1:
                dy = 1

            # calculate the angle of the possible movement
            angle = np.arctan2(dy, dx)

            # calculate the angular difference
            angle_difference = np.arctan2(
                np.sin(angle - previous_angle),
                np.cos(angle - previous_angle)
            )

            # calculate the Gaussian memory weight
            memory_weight = np.exp(
                -(angle_difference ** 2) /
                (2 * self.params["n_sigma"] ** 2)
            )

            memory_weights.append(memory_weight)

        # convert weights into a numpy array
        memory_weights = np.array(memory_weights)

        # normalize weights so probabilities sum to one
        memory_probabilities = memory_weights / np.sum(memory_weights)

        return neighbors, memory_probabilities

    def calculate_movement_probabilities(self):

        # calculate probabilities from the chemical gradient
        neighbors, gradient_probabilities = (
            self.calculate_gradient_probabilities()
        )

        # calculate probabilities from directional memory
        _, memory_probabilities = (
            self.calculate_memory_probabilities()
        )

        # combine the two probability distributions
        combined_weights = gradient_probabilities * memory_probabilities

        # normalize the combined weights
        movement_probabilities = (
            combined_weights / np.sum(combined_weights)
        )

        return neighbors, movement_probabilities

    def move_chemotactically(self):

        # calculate possible movements and their probabilities
        neighbors, probabilities = self.calculate_movement_probabilities()

        # randomly select a movement based on its probability
        selected_index = np.random.choice(
            len(neighbors),
            p=probabilities
        )

        # retrieve the selected neighboring position
        new_x, new_y = neighbors[selected_index]

        # retrieve the current position
        old_x, old_y = self.grid.get_location(self)

        # calculate the actual movement direction
        dx = new_x - old_x
        dy = new_y - old_y

        # correct for periodic boundaries
        if dx > 1:
            dx = -1
        elif dx < -1:
            dx = 1

        if dy > 1:
            dy = -1
        elif dy < -1:
            dy = 1

        # identify the selected direction
        self.prev_move = DIRECTIONS.index((dx, dy))

        # update the neutrophil's position
        self.grid.move_agent(self, new_x, new_y)

    def release_pro(self):
        """
        probabilistically release pro-inflammatory mediator
        at the neutrophil's current location.
        """

        # generate a random number between 0 and 1
        release_probability = random.random()

        # check whether the neutrophil releases mediator
        if release_probability < self.params["prob_neutroReleasePro"]:

            # retrieve the neutrophil's current location
            x, y = self.grid.get_location(self)

            # retrieve the amount of mediator to release
            amount = self.params["inc_neutroReleasePro"]

            # increase the local pro-inflammatory concentration
            add_pro_inflammatory(
                self.environment[0],
                x,
                y,
                amount
            )

    def update_lifespan(self):
        """
        decrease the neutrophil's remaining lifespan by one tick.
        """

        # decrease the remaining lifespan
        self.lifespan -= 1

        # prevent the lifespan from becoming negative
        if self.lifespan < 0:
            self.lifespan = 0

    def become_apoptotic(self):
        """
        replace an active neutrophil with an apoptotic neutrophil
        at the same grid location.
        """

        # retrieve the neutrophil's current location
        x, y = self.grid.get_location(self)

        # create a new apoptotic neutrophil at the same location
        apoptotic_neutrophil = ApoptoticNeutrophil(x, y)

        # remove the active neutrophil from the grid
        self.grid.remove_agent(self)

        # add the apoptotic neutrophil to the grid
        self.grid.add_agent(apoptotic_neutrophil, x, y)

        # decrease the active neutrophil count
        Neutrophil.neutrophil_counter -= 1

        # return the new apoptotic neutrophil
        return apoptotic_neutrophil

    def run(self):

        # check whether the neutrophil is still alive
        if self.lifespan > 0:

            # allow the neutrophil to move multiple times per tick
            for _ in range(self.params["movesPerTick_n"]):
                self.move_chemotactically()

            # release pro-inflammatory mediator
            self.release_pro()

            # decrease the remaining lifespan
            self.update_lifespan()

        else:

            # convert the neutrophil into an apoptotic cell
            self.become_apoptotic()