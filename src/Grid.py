class Grid:

    def __init__(self, width, height):
        """
        initialize the grid and create the data structures
        needed to keep track of all agents.
        """

        # store the width and height of the grid
        self.width = width
        self.height = height

        # create a dictionary that stores all agents at each grid location
        # the key is the (x, y) coordinate
        # the value is a list containing all agents at that location
        # multiple agents are allowed to occupy the same location
        self.cells = {}

        # create a dictionary that stores the location of each agent
        # the key is the agent object
        # the value is the agent's (x, y) coordinate
        # this allows us to quickly find an agent's location
        self.agent_locations = {}


    def wrap_coordinates(self, x, y):
        """
        apply periodic boundary conditions to the grid coordinates.

        if an agent moves beyond one edge of the grid, it reappears
        on the opposite edge.
        """

        # use the modulo operator to keep x within the grid boundaries
        x = x % self.width

        # use the modulo operator to keep y within the grid boundaries
        y = y % self.height

        # return the updated coordinates
        return x, y


    def add_agent(self, agent, x, y):
        """
        add an agent to the grid at the specified coordinates.
        """

        # make sure the coordinates are within the grid boundaries
        x, y = self.wrap_coordinates(x, y)

        # check whether the agent is already present on the grid
        # this prevents the same agent from being added more than once
        if agent in self.agent_locations:
            raise ValueError("agent is already on the grid.")

        # check whether this grid location already exists in the dictionary
        # if it does not exist, create an empty list for this location
        if (x, y) not in self.cells:
            self.cells[(x, y)] = []

        # add the agent to the list of agents at this location
        self.cells[(x, y)].append(agent)

        # record the agent's location in the location dictionary
        self.agent_locations[agent] = (x, y)

        # update the agent's own coordinates
        agent.x = x
        agent.y = y


    def remove_agent(self, agent):
        """
        remove an agent from the grid.
        """

        # check whether the agent is currently present on the grid
        # if it is not present, there is nothing to remove
        if agent not in self.agent_locations:
            return

        # retrieve the agent's current coordinates
        x, y = self.agent_locations[agent]

        # remove the agent from the list at its current location
        self.cells[(x, y)].remove(agent)

        # check whether any agents remain at this location
        # if the list is empty, remove the location from the dictionary
        if len(self.cells[(x, y)]) == 0:
            del self.cells[(x, y)]

        # remove the agent from the location dictionary
        del self.agent_locations[agent]


    def move_agent(self, agent, new_x, new_y):
        """
        move an agent from its current location to a new location.
        """

        # check whether the agent is currently present on the grid
        if agent not in self.agent_locations:
            raise ValueError("cannot move an agent that is not on the grid.")

        # apply periodic boundary conditions to the new coordinates
        new_x, new_y = self.wrap_coordinates(new_x, new_y)

        # remove the agent from its previous grid location
        self.remove_agent(agent)

        # add the agent to its new grid location
        self.add_agent(agent, new_x, new_y)


    def get_location(self, agent):
        """
        return the current coordinates of a specific agent.
        """

        # check whether the agent is present on the grid
        if agent not in self.agent_locations:
            raise ValueError("agent is not on the grid.")

        # return the coordinates stored for this agent
        return self.agent_locations[agent]


    def get_agents_at(self, x, y):
        """
        return all agents occupying a specific grid location.
        """

        # make sure the coordinates are within the grid boundaries
        x, y = self.wrap_coordinates(x, y)

        # retrieve the list of agents at this location
        # if the location has no agents, return an empty list
        return self.cells.get((x, y), [])


    def get_neighbors(self, x, y):
        """
        return the eight neighboring locations surrounding a given location.

        these eight locations form the moore neighborhood.
        """

        # define the eight neighboring directions in a consistent order
        directions = [
            (1, 0),
            (1, 1),
            (0, 1),
            (-1, 1),
            (-1, 0),
            (-1, -1),
            (0, -1),
            (1, -1)
        ]

        # create an empty list to store the neighboring coordinates
        neighbors = []

        # examine each of the eight directions
        for dx, dy in directions:

            # calculate the coordinates of the neighboring location
            neighbor_x = x + dx
            neighbor_y = y + dy

            # apply periodic boundary conditions
            neighbor_x, neighbor_y = self.wrap_coordinates(
                neighbor_x, neighbor_y
            )

            # add the neighboring coordinates to the list
            neighbors.append((neighbor_x, neighbor_y))

        # return the list containing all eight neighboring locations
        return neighbors

    def get_neighborhood_agents(self, x, y):
        """
        return all agents located in the eight neighboring grid locations.

        the agents at the current location are not included.
        """

        # create an empty list to store neighboring agents
        neighborhood_agents = []

        # retrieve the eight neighboring grid locations
        neighbors = self.get_neighbors(x, y)

        # examine each neighboring location
        for neighbor_x, neighbor_y in neighbors:

            # retrieve all agents occupying the current neighboring location
            agents = self.get_agents_at(neighbor_x, neighbor_y)

            # add these agents to the neighborhood list
            # extend adds each agent individually rather than adding a nested list
            neighborhood_agents.extend(agents)

        # return all agents found in the neighboring locations
        return neighborhood_agents


    def get_all_agents(self):
        """
        return a list containing every agent currently on the grid.
        """

        # retrieve all agent objects stored as keys in the location dictionary
        # convert the dictionary keys into a list
        return list(self.agent_locations.keys())


    def get_agents_by_type(self, agent_type):
        """
        return all agents belonging to a specified agent class.

        for example, this can be used to retrieve all neutrophils
        or all macrophages currently on the grid.
        """

        # create an empty list to store agents of the requested type
        matching_agents = []

        # retrieve every agent currently on the grid
        for agent in self.get_all_agents():

            # check whether the agent belongs to the requested class
            if isinstance(agent, agent_type):

                # add the matching agent to the list
                matching_agents.append(agent)

        # return all agents belonging to the requested class
        return matching_agents


    def get_random_neighbor(self, x, y):
        """
        randomly select and return one of the eight neighboring locations.
        """

        # import the random module to select a neighbor
        import random

        # retrieve the eight neighboring locations
        neighbors = self.get_neighbors(x, y)

        # randomly select one neighboring location and return it
        return random.choice(neighbors)


