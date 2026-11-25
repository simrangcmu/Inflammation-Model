#non-chronic inflammation: apoptotic neutrophils accumulate locally (concentrated at site)
    #- recruited macrophages, engulf many - reaches threshold and so back to normal state
#chronic inflammation: macrophages fail to engulf enough apoptotic neutrophils
#- never switches to anti-inflammatory state and the signal becomes too weak so keeps running

class ApoptoticNeutrophil: 
    # track the total number of apoptotic neutrophils
    apoptotic_counter = 0

    def __init__(self, x, y):
        
        #initial state of the apoptotic neutrophil
        self.x = x
        self.y = y
        self.state = "apoptotic"
        self.engulfed_by_macrophage = False
        self.threshold_reached = False
        # increment the total apoptotic neutrophil count
        ApoptoticNeutrophil.apoptotic_counter += 1


    #methods
    def detect_engulfment(self):
        #check to see if available to engulf by a macrophage
        pass

    def is_apoptotic(self):
        """
        check whether the neutrophil is apoptotic.
        """
        return self.state == "apoptotic"

    def is_engulfed(self):
        """
        check whether a macrophage has engulfed the neutrophil.
        """
        return self.engulfed_by_macrophage

    def get_engulfed(self):
        """
        mark the neutrophil as engulfed.
        """
        self.engulfed_by_macrophage = True

    def update_state(self): 
        #code to update the state of the apoptotic neutrophil based on whether it has been engulfed by a macrophage
        pass

#no move method because they stay in the same spot & don't secrete signals  