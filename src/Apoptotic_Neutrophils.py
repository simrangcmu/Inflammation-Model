#non-chronic inflammation: apoptotic neutrophils accumulate locally (concentrated at site)
    #- recruited macrophages, engulf many - reaches threshold and so back to normal state

#chronic inflammation: macrophages fail to engulf enough apoptotic neutrophils
#- never switches to anti-inflammatory state and the signal becomes too weak so keeps running


#the paper's model does not require macrophages to reach a phagocytosis threshold. A macrophage switches to its
#  anti-inflammatory phenotype after its first successful phagocytosis. Subsequent phagocytosis events can continue clearing apoptotic neutrophils.
import numpy as np

from environment import add_pro_inflammatory


class ApoptoticNeutrophil:

    # track the number of apoptotic neutrophils currently present
    apoptotic_counter = 0

    def __init__(self, x, y):

        # assign a random lifespan between 60 and 720 ticks
        self.lifespan = np.random.randint(60, 721)

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

    def undergo_necrosis(self, grid, c, inc_apoptoticReleasePro):
        """
        release pro-inflammatory mediator when the apoptotic
        neutrophil undergoes necrosis, then remove it from the grid.
        """

        # increase the pro-inflammatory concentration at the current location
        add_pro_inflammatory(
            c,
            self.x,
            self.y,
            inc_apoptoticReleasePro
        )

        # remove the apoptotic neutrophil from the grid
        grid.remove_agent(self)

        # update the apoptotic neutrophil count
        ApoptoticNeutrophil.apoptotic_counter -= 1

        # update the state
        self.state = "necrotic"

    def update(self, grid, c, inc_apoptoticReleasePro):
        """
        update the apoptotic neutrophil's lifespan.

        if its lifespan expires, the neutrophil undergoes necrosis.
        """

        # do not update a neutrophil that has already been engulfed
        if self.is_engulfed():
            return

        # check whether the lifespan has expired
        if self.lifespan <= 0:

            # undergo necrosis and remove the neutrophil
            self.undergo_necrosis(
                grid,
                c,
                inc_apoptoticReleasePro
            )

            return

        # decrease the remaining lifespan
        self.lifespan -= 1

#no move method because they stay in the same spot & don't secrete signals  