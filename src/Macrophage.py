#non chronic resolution: macrophages are recruited to them and phagocytose them and 
#become anti-inflammatory

#chronic inflammation: macrophages don't encounter and phagocytose enough apoptoticneutrophils to become activated
#there isn't enough anti-inflammatory signaling to resolve the resposne

"""
Python translation of the Repast Simphony Macrophage agent
(Bayani, Dunster, Crofts, Nelson - inflammatory response model).
 
"""
import math
import random
import numpy as np

from apoptotic_neutrophils import ApoptoticNeutrophil
from environment import add_anti_inflammatory
 
# The 8 Moore-neighbourhood moves: 0=right, 1=up/right, 2=up, 3=up/left, 4=left, 5=down/left, 6=down, 7=down/right


macrophage_params = {
    "movesPerTick_m": 2,

    "prob_macroPhagoApop": 1.0,
    "threshold_macroActivation": 2,

    "prob_macroReleaseAnti": 1.0,
    "inc_macroReleaseAnti": 0.5,

    "m_k_grad": 0.5,

    "prob_macroLeave": 0.0,

    "macroRecruitmentRadius": 10,
}

class Macrophage:

    # track the number of macrophages
    macrophage_counter = 0

    def __init__(self, grid, environment, params):

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

    def find_nearest_apoptotic(self):
        """
        Find the nearest unengulfed apoptotic neutrophil
        within the macrophage recruitment radius.
        """

        x, y = self.grid.get_location(self)

        recruitment_radius = self.params["macroRecruitmentRadius"]

        nearest_apoptotic = None
        shortest_distance = None

        apoptotic_neutrophils = self.grid.get_agents_by_type(
            ApoptoticNeutrophil
        )

        for apoptotic_neutrophil in apoptotic_neutrophils:

            if apoptotic_neutrophil.is_engulfed():
                continue

            apoptotic_x, apoptotic_y = self.grid.get_location(
                apoptotic_neutrophil
            )

            # calculate the distance in each direction
            distance_x = abs(apoptotic_x - x)
            distance_y = abs(apoptotic_y - y)

            # account for periodic boundary conditions
            distance_x = min(distance_x, self.grid.width - distance_x)
            distance_y = min(distance_y, self.grid.height - distance_y)

            # calculate the shortest distance to the apoptotic neutrophil
            distance = math.sqrt(
                distance_x ** 2 + distance_y ** 2
)

            if distance <= recruitment_radius:

                if shortest_distance is None or distance < shortest_distance:

                    shortest_distance = distance
                    nearest_apoptotic = apoptotic_neutrophil

        return nearest_apoptotic


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
        """Change the macrophage activation state after its first successful phagocytosis"""

        # check whether the macrophage has reached its activation threshold
        if self.phagocytosis_count >= self.params["threshold_macroActivation"]:

            # switch the macrophage to its activated state
            self.activated = True

        #to switch it to activated state to relesase anti inflamm mediator

        #WE TEMPORARARILY USE "prob_macroPhagoApop": 1.0, and "threshold_macroActivation": 2
        #CHANGE LATER IF NEEDED


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


    def calculate_movement_probabilities(self):
        """
        Determine movement probabilities based on the pro-inflammatory
        concentration at neighboring locations.
        """

        # get the macrophage's current coordinates
        x, y = self.grid.get_location(self)

        # get the current pro-inflammatory concentration
        current_c = self.environment[0][y, x]

        # get the eight neighboring locations
        neighboring_locations = self.grid.get_neighbors(x, y)

        # create a list to store movement exponents
        movement_exponents = []

        # calculate an exponent for each neighboring location
        for neighbor_x, neighbor_y in neighboring_locations:

            # get the pro-inflammatory concentration at this location
            neighbor_c = self.environment[0][neighbor_y, neighbor_x]

            # calculate how much higher or lower the concentration is
            concentration_difference = neighbor_c - current_c

            # calculate the movement exponent
            exponent = self.params["m_k_grad"] * concentration_difference

            movement_exponents.append(exponent)

        # subtract the largest exponent to improve numerical stability
        maximum_exponent = max(movement_exponents)

        # convert the exponents into movement weights
        movement_weights = [
            np.exp(exponent - maximum_exponent)
            for exponent in movement_exponents
        ]

        # calculate the total movement weight
        total_weight = sum(movement_weights)

        # normalize the weights into probabilities
        movement_probabilities = [
            weight / total_weight
            for weight in movement_weights
        ]

        return movement_probabilities

        #when a macrophage cannot find apoptotic neutrophils nearby, it should 
        # move toward regions with stronger pro-inflammatory signaling.
    
    def move_chemotactically(self):
        """
        Move the macrophage toward higher concentrations
        of pro-inflammatory mediator c.
        """

        # get the macrophage's current coordinates
        x, y = self.grid.get_location(self)

        # get the neighboring locations
        neighboring_locations = self.grid.get_neighbors(x, y)

        # calculate the probability of moving to each location
        movement_probabilities = self.calculate_movement_probabilities()

        # randomly select a neighboring location using the probabilities
        selected_index = np.random.choice(
            len(neighboring_locations),
            p=movement_probabilities
        )

        # get the selected neighboring coordinates
        new_x, new_y = neighboring_locations[selected_index]

        # move the macrophage
        self.grid.move_agent(self, new_x, new_y)

    def move(self):
        """
        Move toward apoptotic neutrophils when they are
        within the macrophage recruitment radius.

        If no apoptotic neutrophils are detected,
        move chemotactically toward higher concentrations
        of pro-inflammatory mediator c.
        """

        # find the nearest apoptotic neutrophil within
        # the macrophage's recruitment radius
        target = self.find_nearest_apoptotic()

        # if no apoptotic neutrophil is nearby,
        # follow the pro-inflammatory gradient
        if target is None:
            self.move_chemotactically()
            return

        # get the macrophage's current location
        macrophage_x, macrophage_y = self.grid.get_location(self)

        # get the apoptotic neutrophil's location
        target_x, target_y = self.grid.get_location(target)

        # calculate the direct distance in each direction
        difference_x = target_x - macrophage_x
        difference_y = target_y - macrophage_y

        # account for periodic boundaries in the x direction
        if abs(difference_x) > self.grid.width / 2:

            if difference_x > 0:
                difference_x -= self.grid.width
            else:
                difference_x += self.grid.width

        # account for periodic boundaries in the y direction
        if abs(difference_y) > self.grid.height / 2:

            if difference_y > 0:
                difference_y -= self.grid.height
            else:
                difference_y += self.grid.height

        # determine the direction of movement
        if difference_x > 0:
            movement_x = 1
        elif difference_x < 0:
            movement_x = -1
        else:
            movement_x = 0

        if difference_y > 0:
            movement_y = 1
        elif difference_y < 0:
            movement_y = -1
        else:
            movement_y = 0

        # calculate the macrophage's next position
        new_x = macrophage_x + movement_x
        new_y = macrophage_y + movement_y

        # move the macrophage one grid cell toward the target
        self.grid.move_agent(self, new_x, new_y)

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