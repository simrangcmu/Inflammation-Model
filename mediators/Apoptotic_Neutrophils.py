#non-chronic inflammation: apoptotic neutrophils accumulate locally (concentrated at site)
    #- recruited macrophages, engulf many - reaches threshold and so back to normal state
#chronic inflammation: macrophages fail to engulf enough apoptotic neutrophils
#- never switches to anti-inflammatory state and the signal becomes too weak so keeps running

class ApoptoticNeutrophil: 
    def __init__(self, x, y):
        #initial state of the apoptotic neutrophil
        self.x = x
        self.y = y
        self.state = "apoptotic"
        self.engulfed_by_macrophage = False
        self.threshold_reached = False


    #methods
    def detect_engulfment(self):
        #check to see if available to engulf by a macrophage
        pass

    def get_engulfed(self):
        #code to mark the apoptotic neutrophil as engulfed by a macrophage
        self.engulfed_by_macrophage = True
        pass

    def is_engulfed(self):
        #code to check if the apoptotic neutrophil has been engulfed by a macrophage
        pass

    def is_apoptotic(self):
        #code to check if the neutrophil is apoptotic
        pass

    def update_state(self): 
        #code to update the state of the apoptotic neutrophil based on whether it has been engulfed by a macrophage
        pass

#no move method because they stay in the same spot & don't secrete signals  