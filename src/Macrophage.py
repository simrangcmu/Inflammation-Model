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
from collections import defaultdict

from apoptotic_neutrophils import ApoptoticNeutrophil
from environment import add_anti_inflammatory
 
# The 8 Moore-neighbourhood moves: 0=right, 1=up/right, 2=up, 3=up/left, 4=left, 5=down/left, 6=down, 7=down/right

Directions = [(0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (1, 0), (1, 1)]
 
params = {
    "movesPerTick_n": 2,
    "n_k_grad": 0.5,
    "n_sigma": 1.0,
    "prob_neutroReleasePro": 0.5,
    "inc_neutroReleasePro": 0.1,
}

macrophage_params = {
    "movesPerTick_m": 2,

    "prob_macroPhagoApop": 1.0,
    "threshold_macroActivation": 2,

    "prob_macroReleaseAnti": 1.0,
    "inc_macroReleaseAnti": 0.5,

    "m_k_grad": 0.5,

    "thresh_macroLeave": 0.5,
    "prob_macroLeave": 0.8
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

        # initialize previous movement direction
        self.prev_move = random.randint(0, 7)

        # macrophages start in an inactive state
        self.activated = False

        # count successful phagocytosis events
        self.phagocytosis_count = 0

        # track whether the macrophage is still alive
        self.alive = True

        # increment the macrophage count
        Macrophage.macrophage_counter += 1

    def find_apoptotic_neighbors(self):
        """
        Identify apoptotic neutrophils in the current location
        and surrounding Moore neighborhood.
        """

        # get the macrophage's current coordinates
        x, y = self.grid.get_location(self)

        # get agents in the current location
        current_agents = self.grid.get_agents_at(x, y)

        # get agents in the eight surrounding locations
        neighboring_agents = self.grid.get_neighborhood_agents(x, y)

        # combine both lists
        nearby_agents = current_agents + neighboring_agents

        # create a list to store apoptotic neutrophils
        apoptotic_neighbors = []

        # check each nearby agent
        for agent in nearby_agents:

            # check whether the agent is an apoptotic neutrophil
            if isinstance(agent, ApoptoticNeutrophil):

                # make sure it has not already been engulfed
                if not agent.is_engulfed():
                    apoptotic_neighbors.append(agent)

        return apoptotic_neighbors


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
    def move_toward_apoptotics(self):
        x, y = self.grid.get_location(self)
 
        # Initialize list of current cell + 8 neighbours 
        candidates = [(x, y)] + self.grid.moore_neighbour_points(x, y)
        random.shuffle(candidates)  # remove bias when counts tie
 
        best_point, max_count = None, -1

        #determine which grid sqaure has most neutrophils

        for point in candidates:
            count = 0
            for obj in self.grid.objects_at(x, y): #need to add grid.objects_at() getter function
                if isinstance(obj, ApoptoticNeutrophil):#isinstance(obj, Apoptotic) returns True if there is an apoptotic there, and False if not
                    count += 1
            if count > max_count:   
                best_point, max_count = point, count

        #if there are no apoptotic cells in neighborhood, cell will move chemotactically towards pro-inflam mediators, 
        # else move towards space with most apoptotics
        if max_count <= 0:
            self.move_chemotactically()
        elif best_point != (x, y):
            self.grid.move_to(self, *best_point) #use of asterisk unpacks tuple best_point into separate args
 
    #Phagocytose apoptotic neutrophils

    def phagocytose_apoptotic(self):
        
        x, y = self.grid.get_location(self)
        #make list of apoptotic agents of current grid square
        apoptotics = []
        for obj in self.grid.objects_at(x, y):
            if isinstance(obj, ApoptoticNeutrophil):
                apoptotics.append(obj)
        
        if len(apoptotics) > 0:
            target = random.choice(apoptotics)  #chooses random object from list, with equal weighted probabilities
            if random.random() < self.params["prob_macroPhagoApop"]:
                target.die()
                self.release_anti_on = True


    #Release of anti-inflammatory mediators
    def release_anti(self):
        if self.release_anti_on:
            x, y = self.grid.get_location(self)
            here = self.grid.environment_at(x, y)
            if random.random() < self.params["prob_macroReleaseAnti"]:
                here.increase_anti_concentration(self.params["inc_macroReleaseAnti"])

    #Make decision of if agents will leave tissue/grid
    def make_leave_decision(self):
        x, y = self.grid.get_location(self)
        local_pro = 0.0
        for pt in self.grid.moore_neighbour_points(x, y):
            env = self.grid.environment_at(*pt)
            if env is not None:
                local_pro += env.pro_concentration
 
        if local_pro < self.params["thresh_macroLeave"]:
            if random.random() < self.params["prob_macroLeave"]:
                self.die()
 
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
        """Determine movement probabilites based on apoptotic-cell attraction 
         or inflammatory conditions
         """

        # get the macrophage's current coordinates
        x, y = self.grid.get_location(self)

        # get the current pro-inflammatory concentration
        current_c = self.environment[0][y, x]

        # get the eight neighboring locations
        neighboring_locations = self.grid.get_neighbors(x, y)

        # create a list to store movement weights
        movement_weights = []

        # calculate a weight for each neighboring location
        for neighbor_x, neighbor_y in neighboring_locations:

            # get the pro-inflammatory concentration at this location
            neighbor_c = self.environment[0][neighbor_y, neighbor_x]

            # calculate how much higher or lower the concentration is
            concentration_difference = neighbor_c - current_c

            # assign a higher weight to locations with higher c
            weight = np.exp(
                self.params["m_k_grad"] * concentration_difference
            )

            movement_weights.append(weight)

        # convert movement weights into probabilities
        total_weight = sum(movement_weights)

        movement_probabilities = [
            weight / total_weight for weight in movement_weights
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
        Move toward nearby apoptotic neutrophils.
        If none are nearby, move chemotactically toward
        higher concentrations of pro-inflammatory mediator c.
        """

        # identify apoptotic neutrophils in the neighborhood
        apoptotic_neighbors = self.find_apoptotic_neighbors()

        # check whether any apoptotic neutrophils were found
        if len(apoptotic_neighbors) == 0:

            # move toward higher concentrations of c
            self.move_chemotactically()

        else:

            # get the macrophage's current coordinates
            x, y = self.grid.get_location(self)

            # get the eight neighboring locations
            neighboring_locations = self.grid.get_neighbors(x, y)

            # include the current location
            candidate_locations = [(x, y)] + neighboring_locations

            # shuffle locations to avoid directional bias
            random.shuffle(candidate_locations)

            # find the location containing the most apoptotic neutrophils
            best_location = None
            highest_count = 0

            for location in candidate_locations:

                candidate_x, candidate_y = location

                agents = self.grid.get_agents_at(candidate_x, candidate_y)

                apoptotic_count = 0

                for agent in agents:

                    if isinstance(agent, ApoptoticNeutrophil):

                        if not agent.is_engulfed():
                            apoptotic_count += 1

                if apoptotic_count > highest_count:

                    highest_count = apoptotic_count
                    best_location = location

            # move toward the selected location
            if best_location is not None and best_location != (x, y):

                new_x, new_y = best_location

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

    def check_death(self):
        """
        Check whether the macrophage has reached the end of its lifespan.
        """

        # check whether the macrophage has expired
        if self.lifespan <= 0:

            # remove the macrophage
            self.die()

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
            self.move()

        # attempt to engulf an apoptotic neutrophil at its current location
        self.phagocytose()

        # release anti-inflammatory mediator if the macrophage is activated
        self.release_mediator()

        # check whether the macrophage should leave the tissue
        self.check_departure()

        # stop if the macrophage has left the tissue
        if not self.alive:
            return

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