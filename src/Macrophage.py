#non chronic resolution: macrophages are recruited to them and phagocytose them and 
#become anti-inflammatory

#chronic inflammation: macrophages don't encounter and phagocytose enough apoptoticneutrophils to become activated
#there isn't enough anti-inflammatory signaling to resolve the resposne

class Macrophage:
    def __init__(self):
        #initial state of the macrophage
        self.state = "resting"
        self.apoptotic_neutrophils_engulfed = 0
        self.signaling_molecules_secreted = 0


    #macrophage parameters
        #add here 

    #methods
    def phagocytose_apoptotic_neutrophils(self):
        #code to phagocytose apoptotic neutrophils 
        pass

    def become_anti_inflammatory(self): 
        #threshold to become anti-inflammatory
        pass

#behaviors of macrophage
    def move(self): 
        #code to move the macrophage to engulf apoptotic neutrophils
        pass 

    def detect_apoptotic_neutrophils(self): 
        #code to detect apoptotic neutrophils
        pass

    def engulf_apoptotic_neutrophils(self):
        #code to engulf apoptotic neutrophils
        pass

    def update_state(self): 
        #code to update the state of the macrophage based on interactions with apoptotic neutrophils
        pass

    def secrete_signaling_molecules(self): 
        #code to secrete signaling molecules based on the state of the macrophage
        pass

    def pro_inflammatory_macrophage(self): 
        #code to represent the pro-inflammatory state of the macrophage
        pass

    def anti_inflammatory_macrophage(self):
        #code to represent the anti-inflammatory state of the macrophage
        pass

    def fixed_type_macrophage(self):
        #code to represent the fixed position of the macrophage
        pass

    def free_type_macrophage(self):
        #code to represent the macrophage that can move
        pass   
