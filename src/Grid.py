
    def wrap_coordinates(self, x: int, y: int):
        """Maintain boundary of grid
        Returns: tuple containing wrapped x and y coordinates"""

        x = x % self.width
        y = y % self.height


    def add_agent():
        """Place agent at cell and update its position
        
        calls wrap_coordinates() to check position?"""

    def remove_agent():
        """ Remove agent from current position"""

    def move_agent(self, agent, new_x, new_y):
        """Remove agent from old position, and move to new position"""

        #Wrap new coordinates around grid
        new_x, new_y = self.wrap_coordinates(new_x, new_y)

        #Remove agent from old position
        self.remove_agent(agent)

        #Add agents at new position
        self.add_agent(agent, new_x, new_y)


    def get_agents_at():
        """Return all agents at given cell"""

    def get_neighbors():
        """Return 8 Moore neighborhood coordinates given center cell"""

    def get_neighborhood_agents():
        """ Return agents occupying 8 Moore neighborhood cells"""

    def get_all_agents():
        """ Return all agents currently on grid"""

    def get_agents_by_type():
        """Return agents of given type (macrophage, neutrophil, apoptotic)"""

    def get_random_neighbor():
        """Return randomly selected neighboring coordinate"""

    
