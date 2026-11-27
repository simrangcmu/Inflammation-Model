#non-chronic inflammation: apoptotic neutrophils accumulate locally (concentrated at site)
    #- recruited macrophages, engulf many - reaches threshold and so back to normal state
#chronic inflammation: macrophages fail to engulf enough apoptotic neutrophils
#- never switches to anti-inflammatory state and the signal becomes too weak so keeps running

class ApoptoticNeutrophil:

    # track the number of apoptotic neutrophils currently present
    apoptotic_counter = 0

    def __init__(self, x, y):

        # initialize the apoptotic neutrophil
        self.x = x
        self.y = y
        self.state = "apoptotic"
        self.engulfed_by_macrophage = False

        # increment the apoptotic neutrophil count
        ApoptoticNeutrophil.apoptotic_counter += 1

    def is_apoptotic(self):
        """
        Check whether the neutrophil is apoptotic.
        """
        return self.state == "apoptotic"

    def is_engulfed(self):
        """
        Check whether a macrophage has engulfed the neutrophil.
        """
        return self.engulfed_by_macrophage

    def get_engulfed(self):
        """
        Mark the neutrophil as engulfed.
        """
        self.engulfed_by_macrophage = True
        self.state = "engulfed"



    #methods
    def detect_engulfment(self):
        #check to see if available to engulf by a macrophage
        pass


    def update_state(self): 
        #code to update the state of the apoptotic neutrophil based on whether it has been engulfed by a macrophage
        pass

#no move method because they stay in the same spot & don't secrete signals  