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
 
# The 8 Moore-neighbourhood moves: 0=right, 1=up/right, 2=up, 3=up/left, 4=left, 5=down/left, 6=down, 7=down/right

Directions = [(0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (1, 0), (1, 1)]
 
params = {"movesPerTick_n": 2, "n_k_grad": 0.5, "n_sigma": 1.0,
    "prob_neutroReleasePro": 0.5, "inc_neutroReleasePro": 0.1,}

class Macrophage:
    macrophage_counter = 0  #initialize total count of macrophages to 0 
 
    def __init__(self, grid, params):
        self.grid = grid
        self.params = params  
        self.lifespan = np.random.Generator(1440, 86400)   # initializes lifespan; using numpy because of recommendation of better reproducibility than random.randint 
        self.prev_move = np.random.Generator(0, 7) #initializes random previous move
        self.release_anti_on = False #marcrophage does not release until it has phagocytosed an apoptotic agent
        Macrophage.macrophage_counter += 1 #increments macrophage_counter eery time new Marcrophage object is created
 

    def run(self):
        if self.lifespan > 0:
            for _ in range(self.params["movesPerTick_m"]): #what is the "movesPerTick_m" parameter mean??
                self.move_toward_apoptotics()
            
            #try to phagocytose apoptotic neutrophils
            self.phagocytose_apoptotic()
            
            #release anti-inflam mediator if conditions met
            self.release_anti()

            #decide if vacate tissue
            self.make_leave_decision()

            #decrement lifespan
            self.lifespan -= 1
        else: #if lifepsan = 0, the cell dies
            self.die()
 
    def die(self):
        Macrophage.macrophage_counter -= 1
        self.grid.remove(self) #attribute of grid class
 

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
                if isinstance(obj, Apoptotic):#isinstance(obj, Apoptotic) returns True if there is an apoptotic there, and False if not
                    count += 1
            if count > max_count:   
                best_point, max_count = point, count

        #if there are no apoptotic cells in neighborhood, cell will move chemotactically towards pro-inflam mediators, 
        # else move towards space with most apoptotics
        if max_count <= 0:
            self.move_chemotactically()
        elif best_point != (x, y):
            self.grid.move_to(self, *best_point) #use of asterisk unpacks tuple best_point into separate args
 

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
 
    #Phagocytose apoptotic neutrophils

    def phagocytose_apoptotic(self):
        
        x, y = self.grid.get_location(self)
        #make list of apoptotic agents of current grid square
        apoptotics = []
        for obj in self.grid.objects_at(x, y):
            if isinstance(obj, Apoptotic):
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
 


  
