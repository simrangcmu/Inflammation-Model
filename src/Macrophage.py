#non chronic resolution: macrophages are recruited to them and phagocytose them and 
#become anti-inflammatory

#chronic inflammation: macrophages don't encounter and phagocytose enough apoptoticneutrophils to become activated
#there isn't enough anti-inflammatory signaling to resolve the resposne

import math
import random
import numpy as np

from apoptotic_neutrophils import ApoptoticNeutrophil
from environment import add_anti_inflammatory
 
# The 8 Moore-neighbourhood moves: 0=right, 1=up/right, 2=up, 3=up/left, 4=left, 5=down/left, 6=down, 7=down/right

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

macrophage_params = {
    "movesPerTick_m": 1,

    "prob_macroPhagoApop": 1.0,

    "prob_macroReleaseAnti": 0.8,
    "inc_macroReleaseAnti": 0.001,

    "m_k_grad": 80,
    "m_sigma": 1.0,

    "prob_macroLeave": 0.8,
    "thresh_macroLeave": 0.02,
}

class Macrophage:

    # track the number of macrophages
    macrophage_counter = 0

    def __init__(self, grid, environment, params):

        self.prev_move = None
        self.grid = grid
        self.environment = environment
        self.params = params

        # initialize macrophage lifespan
        self.lifespan = np.random.randint(1440, 86401)

        # macrophages start in an inactive state
        self.activated = False

        # count successful phagocytosis events
        self.phagocytosis_count = 0

        # track whether the macrophage is still alive
        self.alive = True

        # increment the macrophage count
        Macrophage.macrophage_counter += 1


    def die(self):
        """
        Remove the macrophage from the grid when it dies or leaves.
        """

        # make sure the macrophage has not already been removed
        if not self.alive:
            return

        # remove the macrophage from the grid
        self.grid.remove_agent(self)

        # update the macrophage's status
        self.alive = False

        # decrease the macrophage population count
        Macrophage.macrophage_counter -= 1
 

    #governs movement towards apoptotic neutrophils in neighborhood, if any are there


    def find_apoptotic_neighbors(self):
        """
        find apoptotic neutrophils in the eight neighboring locations.
        """

        # get the macrophage's current location
        x, y = self.grid.get_location(self)

        # retrieve the eight neighboring locations
        neighboring_locations = self.grid.get_neighbors(x, y)

        # store locations containing apoptotic neutrophils
        apoptotic_locations = []

        # examine each neighboring location
        for neighbor_x, neighbor_y in neighboring_locations:

            # retrieve all agents at this location
            agents = self.grid.get_agents_at(neighbor_x, neighbor_y)

            # check whether an apoptotic neutrophil is present
            for agent in agents:

                if isinstance(agent, ApoptoticNeutrophil):

                    if agent.is_apoptotic():

                        # store this location only once
                        apoptotic_locations.append((neighbor_x, neighbor_y))
                        break

        return apoptotic_locations

    def find_apoptotic_at_current_location(self):
        """Identity apoptotic cells occupying the macrophages current coordinates"""

        # get the macrophage's current coordinates
        x, y = self.grid.get_location(self)

        # get all agents occupying the same location
        current_agents = self.grid.get_agents_at(x, y)

        # create a list to store apoptotic neutrophils
        apoptotic_cells = []

        # check each agent at the current location
        for agent in current_agents:

            # check whether the agent is an apoptotic neutrophil
            if isinstance(agent, ApoptoticNeutrophil):

                # make sure it has not already been engulfed
                if not agent.is_engulfed():
                    apoptotic_cells.append(agent)

        return apoptotic_cells

    #returns empty list if apoptotic neutrophil is not occupying the macrophages lockation
    #we need it to occupy same location for phagocytosis
    
    def phagocytose(self):
        """Attempt to engulf an apoptotic cell with probability (p_ma);
        remove if it was successful"""

        # find apoptotic neutrophils at the current location
        apoptotic_cells = self.find_apoptotic_at_current_location()

        # check whether there are any apoptotic neutrophils present
        if len(apoptotic_cells) == 0:
            return

        # randomly select one apoptotic neutrophil
        target = random.choice(apoptotic_cells)

        # generate a random number to determine whether phagocytosis succeeds
        phagocytosis_probability = random.random()

        # check whether phagocytosis is successful
        if phagocytosis_probability < self.params["prob_macroPhagoApop"]:

            # mark the apoptotic neutrophil as engulfed
            target.get_engulfed()

            # remove the apoptotic neutrophil from the grid
            self.grid.remove_agent(target)

            # decrease the number of apoptotic neutrophils present
            ApoptoticNeutrophil.apoptotic_counter -= 1

            # increase the macrophage's successful phagocytosis count
            self.phagocytosis_count += 1

            # check whether the macrophage should become activated
            self.activate()

        #if prob is higher than parameter of probabiloty then phagocytosis fails
        #once macrophage engulfs apoptotic neutrophil, we remove it so multiple macrophages dont engulf same cell

    def activate(self):
        """
        activate the macrophage after its first successful phagocytosis.
        """

        # switch to the anti-inflammatory phenotype after the first engulfment
        if self.phagocytosis_count >= 1:
            self.activated = True

    def release_mediator(self):
        """If activated, probabilistically release anti-inflammatory mediator g"""

         # check whether the macrophage is activated
        if not self.activated:
            return

        # generate a random number to determine whether mediator is released
        release_probability = random.random()

        # check whether mediator release is successful
        if release_probability < self.params["prob_macroReleaseAnti"]:

            # get the macrophage's current coordinates
            x, y = self.grid.get_location(self)

            # get the amount of anti-inflammatory mediator to release
            amount = self.params["inc_macroReleaseAnti"]

            # increase the anti-inflammatory concentration at this location
            add_anti_inflammatory(self.environment[1], x, y, amount)


    def calculate_gradient_probabilities(self):
        """
        calculate movement probabilities based on differences
        in pro-inflammatory mediator concentration.
        """

        # retrieve the macrophage's current location
        x, y = self.grid.get_location(self)

        # retrieve the eight neighboring locations
        neighbors = self.grid.get_neighbors(x, y)

        # retrieve the current pro-inflammatory concentration
        current_c = self.environment[0][y, x]

        # initialize a list to store movement weights
        gradient_weights = []

        # calculate the weight for each neighboring location
        for neighbor_x, neighbor_y in neighbors:

            # retrieve the neighboring concentration
            neighbor_c = self.environment[0][neighbor_y, neighbor_x]

            # calculate the concentration difference
            difference = neighbor_c - current_c

            # calculate the gradient weight
            weight = np.exp(
                self.params["m_k_grad"] * difference
                - self.params["m_k_grad"] * max(
                    neighbor_c - current_c
                    for neighbor_x, neighbor_y in neighbors
                    for neighbor_c in [self.environment[0][neighbor_y, neighbor_x]]
                )
)
            gradient_weights.append(weight)

        # convert the weights into a numpy array
        gradient_weights = np.array(gradient_weights)

        # normalize the weights to obtain probabilities
        gradient_probabilities = (
            gradient_weights / np.sum(gradient_weights)
        )

        return neighbors, gradient_probabilities

        #when a macrophage cannot find apoptotic neutrophils nearby, it should 
        # move toward regions with stronger pro-inflammatory signaling.
    
    def calculate_memory_probabilities(self):
        """
        calculate movement probabilities based on the macrophage's
        previous movement direction.
        """

        # retrieve the standard deviation for directional persistence
        sigma = self.params["m_sigma"]

        # if the macrophage has not moved before, assign equal probabilities
        if self.prev_move is None:
            return np.ones(8) / 8

        # convert the previous movement direction into an angle
        previous_angle = self.prev_move * (math.pi / 4)

        # initialize a list to store memory weights
        memory_weights = []

        # calculate the weight for each possible movement direction
        for direction_index in range(8):

            # convert the possible movement direction into an angle
            current_angle = direction_index * (math.pi / 4)

            # calculate the angular difference
            angle_difference = abs(current_angle - previous_angle)

            # account for the circular nature of angles
            angle_difference = min(
                angle_difference,
                2 * math.pi - angle_difference
            )

            # calculate the Gaussian memory weight
            weight = np.exp(
                -(angle_difference ** 2) / (2 * sigma ** 2)
            )

            memory_weights.append(weight)

        # convert the weights into a numpy array
        memory_weights = np.array(memory_weights)

        # normalize the weights to obtain probabilities
        memory_probabilities = (
            memory_weights / np.sum(memory_weights)
        )

        return memory_probabilities


    def calculate_movement_probabilities(self):
        """
        combine the chemical gradient and directional memory
        to calculate the final movement probabilities.
        """

        # calculate the gradient-based movement probabilities
        neighbors, gradient_probabilities = (
            self.calculate_gradient_probabilities()
        )

        # calculate the memory-based movement probabilities
        memory_probabilities = self.calculate_memory_probabilities()

        # combine the gradient and memory probabilities
        combined_weights = (
            gradient_probabilities * memory_probabilities
        )

        # normalize the combined weights
        movement_probabilities = (
            combined_weights / np.sum(combined_weights)
        )

        return neighbors, movement_probabilities
    
    def move_chemotactically(self):
        """
        move the macrophage according to the chemical gradient
        and its previous movement direction.
        """

        # calculate possible movements and their probabilities
        neighbors, probabilities = self.calculate_movement_probabilities()

        # randomly select a neighboring location
        selected_index = np.random.choice(
            len(neighbors),
            p=probabilities
        )

        # retrieve the selected location
        new_x, new_y = neighbors[selected_index]

        # retrieve the current location
        old_x, old_y = self.grid.get_location(self)

        # calculate the movement direction
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

        # remember the direction of this movement
        self.prev_move = DIRECTIONS.index((dx, dy))

        # update the macrophage's location
        self.grid.move_agent(self, new_x, new_y)
     
def move(self):
    """
    Move the macrophage toward apoptotic neutrophils or follow
    the pro-inflammatory mediator gradient.
    """

    # Check whether an apoptotic neutrophil is already at macrophage's current location
    # If so, stay in place so phagocytosis can occur
    if self.find_apoptotic_at_current_location():
        return

    # Find neighboring locations containing apoptotic neutrophils
    apoptotic_neighbors = self.find_apoptotic_neighbors()

    # If apoptotic neutrophils are nearby, move toward one of them
    if apoptotic_neighbors:
        target_x, target_y = random.choice(apoptotic_neighbors)

        dx = target_x - self.grid.get_location(self)[0]
        dy = target_y - self.grid.get_location(self)[1]

        # Account for wrapped boundaries
        if dx > self.grid.width / 2:
            dx -= self.grid.width
        elif dx < -self.grid.width / 2:
            dx += self.grid.width

        if dy > self.grid.height / 2:
            dy -= self.grid.height
        elif dy < -self.grid.height / 2:
            dy += self.grid.height

        self.prev_move = (dx, dy)
        self.grid.move_agent(self, (target_x, target_y))

    # Otherwise, check whether pro-inflammatory mediators are present in neighboring locations
    else:
        x, y = self.grid.get_location(self)
        neighbors = self.grid.get_neighbors(x, y)

        total_c = sum(
            self.environment[0][neighbor_y, neighbor_x]
            for neighbor_x, neighbor_y in neighbors
        )

        # If no pro-inflammatory mediator is present, move randomly
        if total_c == 0:
            target_x, target_y = random.choice(neighbors)

            self.prev_move = (
                target_x - x,
                target_y - y
            )

            self.grid.move_agent(self, (target_x, target_y))

        # Otherwise, follow move chemotactically 
        else:
            self.move_chemotactically()

 
    def check_departure(self):
        """
        Check the local pro-inflammatory concentration and determine
        whether the macrophage leaves the tissue.
        """

        # get the macrophage's current coordinates
        x, y = self.grid.get_location(self)

        # get the eight neighboring locations
        neighboring_locations = self.grid.get_neighbors(x, y)

        # initialize the total local pro-inflammatory concentration
        local_pro = 0.0

        # calculate the total concentration in neighboring locations
        for neighbor_x, neighbor_y in neighboring_locations:

            # get the pro-inflammatory concentration at this location
            local_pro += self.environment[0][neighbor_y, neighbor_x]

        # compare local concentration against the departure threshold
        if local_pro < self.params["thresh_macroLeave"]:

            # generate a random number to determine whether the macrophage leaves
            departure_probability = random.random()

            # check whether the macrophage leaves
            if departure_probability < self.params["prob_macroLeave"]:
                self.die()

            #There are two parameters:
            # thresh_macroLeave: the concentration threshold below which the macrophage may leave
            # prob_macroLeave: the probability that it leaves when the concentration is below that threshold


    def update_lifespan(self):
        """
        Decrease the macrophage's remaining lifespan by one.
        """

        # decrease the remaining lifespan
        self.lifespan -= 1
    
    def update(self):
        """
        Run one complete simulation tick for the macrophage.
        """

        # check whether the macrophage is still alive
        if not self.alive:
            return

        # check whether the macrophage has reached the end of its lifespan
        if self.lifespan <= 0:
            self.die()
            return

        # move the macrophage the required number of times
               
        for _ in range(self.params["movesPerTick_m"]):

            # move toward an apoptotic neutrophil or follow the c gradient
            self.move()

            # attempt phagocytosis at the new location
            self.phagocytose()

            # stop if the macrophage is no longer alive
            if not self.alive:
                return

        # stop if the macrophage has left the tissue
        if not self.alive:
            return

        # release anti-inflammatory mediator if the macrophage is activated
        self.release_mediator()

        # check whether the macrophage should leave the tissue
        self.check_departure()

        # decrease the macrophage's remaining lifespan
        self.update_lifespan()
        
'''
Each time update() is called:
1. Alive check: If the macrophage has already been removed, it does nothing
2. Lifespan check: If its lifespan has expired, it dies
3. Movement: It moves according to movesPerTick_m
4. Phagocytosis: It attempts to engulf an apoptotic neutrophil at its current location
5. Mediator release: If activated, it may release anti-inflammatory mediator
6. Departure: It checks whether local pro-inflammatory signaling is low enough to leave
7. Aging: Its lifespan decreases by one
'''
