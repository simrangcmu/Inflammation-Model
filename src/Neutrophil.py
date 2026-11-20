"""
Python translation of the Repast Simphony Neutrophil agent
(Bayani, Dunster, Crofts, Nelson - inflammatory response model).
 
"""
import math
import random
import numpy as np
from collections import defaultdict
 
# The 8 Moore-neighbourhood moves: 0=right, 1=up/right, 2=up, 3=up/left, 4=left, 5=down/left, 6=down, 7=down/right

Directions = [(0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (1, 0), (1, 1)]
 
params = {"movesPerTick_n": 2, "n_k_grad": 0.5, "n_sigma": 1.0,
    "prob_neutroReleasePro": 0.5, "inc_neutroReleasePro": 0.1,}

class Neutrophil:
    neutrophil_counter = 0  #initialize total count of neutrophil to 0 
 
    def __init__(self, grid, params):
        self.grid = grid
        self.params = params  
        self.lifespan = np.random.Generator(60, 1440)   # initializes lifespan; using numpy because of recommendation of better reproducibility than random.randint 
        self.prev_move = np.random.Generator(0, 7) #initializes random previous move
        Neutrophil.neutrophil_counter += 1 #increments neutrophil_counter eery time new neutrophil object is created
 

    def run(self):
        if self.lifespan > 0:
            for _ in range(self.params["movesPerTick_m"]): #what is the "movesPerTick_m" parameter mean??
                self.move_chemotactically()

            #release pro-inflam mediator
            self.release_pro()

            #decrement lifespan
            self.lifespan -= 1

        else: #if lifepsan = 0, the neutrophil turns into apoptotic cell
            self.become_apoptotic()

    #Release pro-inflam mediator
    def release_pro(self):
            #get current location
            x, y = self.grid.get_location(self)

            #get environment at current location
            here = self.grid.environment_at(x, y)

            #if random value is less than parameter given, increment pro-inflam mediator by parameter amount
            if random.random() < self.params["prob_neutroReleasePro"]:
                here.increase_pro_concentration(self.params["inc_neutroReleasePro"])

    #Remove neutrophil after death
    def die(self):
        Neutrophil.neutrophil_counter -= 1
        self.grid.remove(self)

    #Replace neutrophil with apoptotic agent on apoptosis
    def become_apoptotic(self):
        x, y = self.grid.get_location(self)
        
        # make new apoptotic cell at same spot
        self.grid.add(Apoptotic(self.grid), x, y)  
        
        #remove the original neutrophil 
        self.die()                                  

   #Translate above functions into code structure below

    def get_neighbors():
        """Retrieve 8 neighboring coordinates"""

    def calculate_gradient_probabilities():
        """Calculate movement based on local differences in c"""

    def calculate_memory_probabilities():
        """Calculate directional preference based on previous movement direction, 
        using a Gaussian distribution
        """
    def calculate_movement_probability():
        """ Combine gradient and movement-persistence
        probabilities and normalize result"""

    def move (grid, environment):
        """Select movement direction using the calculated
        probabilities and update position
        """

    def release_mediatiors(environment):
        """ Probabilistically release pro-inflammatory 
        mediator c at current position
        """

    def update_lifespan():
        """ Increment age or decrement remaining lifespan"""

    def check_apoptosis():
        """Determine whether neutrophil has reached
        the end of its lifespan
        """

    def become_apoptotic():
        """Replace active neutrophil with apoptotic neutrophil
        at same grid space/cell
        """

    def update():
        """ Coordinate one active neutrophil's behavior
        with an apoptotic neutrophil at same position
        """
    

    

