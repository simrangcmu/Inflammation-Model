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

 

    #move chemotactically, taking into acount gradient and directional persistence
    def move_chemotactically(self):

        x, y = self.grid.get_location(self)
        here = self.grid.environment_at(x, y)
        neighbour_pts = self.grid.moore_neighbour_points(x, y)  # index == move number
        neighbours = [self.grid.environment_at(*pt) for pt in neighbour_pts]
 
        # 1) Gradient term: weight = exp(k * (neighbour_conc - here_conc))
        k = self.params["m_k_grad"]
        w_grad = [math.exp(k * (n.pro_concentration - here.pro_concentration)) for n in neighbours]
        w_total = sum(w_grad)
        if w_total < 1e-6:
            p_grad = [0.125] * 8
        else:
            p_grad = [w / w_total for w in w_grad]
 
        # 2) Memory term: Gaussian centred on "keep going the same way".
        # shift = how many 45-degree steps this move is from the previous one.
        shift_to_angle = [0, -math.pi / 4, -math.pi / 2, -3 * math.pi / 4,
                          math.pi, 3 * math.pi / 4, math.pi / 2, math.pi / 4]
        sigma2 = self.params["m_sigma"] ** 2
        p_mem = []
        for move in range(8):
            shift = (self.prev_move - move) % 8      # Python % is never negative
            angle = shift_to_angle[shift]
            p_mem.append(math.exp(-angle ** 2 / (2 * sigma2))
                         / math.sqrt(2 * math.pi * sigma2))
        s = sum(p_mem)
        new_list = []
        for p in p_mem:
            new_list.append(p / s)
        p_mem = new_list
 
        # 3) Combine, then pick a move at random using those weights
        # random.choices normalises weights itself, so no manual re-normalising needed
        p_tot = []
        for i in range(8):
            p_tot.append(p_mem[i] * p_grad[i])
        move = random.choices(range(8), weights=p_tot)[0]
 
        self.prev_move = move
        self.grid.move_to(self, *neighbour_pts[move])
 
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


