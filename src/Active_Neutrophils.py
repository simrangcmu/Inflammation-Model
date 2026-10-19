#chronic: neutrophils remain active for too long

class ActiveNeutrophil:
    def __init__(self, x, y):
        #initial state of the active neutrophil
        self.x = x
        self.y = y
        self.state = "active"

    #add parameters for the active neutrophil if needed here 


    #methods
    def is_active(self):
        #code to check if the neutrophil is active
        pass

    def update_state(self):
        #code to update the state of the neutrophil based on its activity
        pass

    def promote_inflammation(self):
        #code to promote inflammation by the active neutrophil
        pass

    def become_apoptotic(self):
        #code to transition the active neutrophil to an apoptotic state
        pass

    def move(self, new_x, new_y):
        #code to move the active neutrophil to a new location --> move towards injury site
        self.x = new_x
        self.y = new_y

    def secrete_signals(self):
        #code for the active neutrophil to secrete inflammatory signals
        pass