import numpy as np
import matplotlib.pyplot as plt


# Grid parameters as given in the paper
GRID_SIZE = 100
DAMAGE_RADIUS = 10
INITIAL_C = 1.0


def initialize_fields():
    """
    initialize the pro-inflammatory and anti-inflammatory mediator fields.

    returns:
        c: pro-inflammatory mediator concentration matrix.
        g: anti-inflammatory mediator concentration matrix.
    """

    # create the two chemical concentration fields
    # every grid location initially has zero concentration
    c = np.zeros((GRID_SIZE, GRID_SIZE), dtype=float)
    g = np.zeros((GRID_SIZE, GRID_SIZE), dtype=float)

    return c, g


def initialize_damage_region(c):
    """
    initialize the circular damaged tissue region.

    the damaged region is assigned an initial pro-inflammatory
    mediator concentration.

    arguments:
        c: pro-inflammatory mediator concentration matrix.

    returns:
        c: updated pro-inflammatory mediator concentration matrix.
    """

    # find the center coordinates of the grid
    center_x = GRID_SIZE // 2
    center_y = GRID_SIZE // 2

    # create coordinates for every grid location
    x, y = np.ogrid[:GRID_SIZE, :GRID_SIZE]

    # calculate the squared distance from the center
    distance_squared = (
        (x - center_x) ** 2 + (y - center_y) ** 2
    )

    # identify the grid locations inside the circular damage region
    damage_region = distance_squared <= DAMAGE_RADIUS ** 2

    # assign the initial pro-inflammatory concentration
    # to every location inside the damaged region
    c[damage_region] = INITIAL_C

    return c


def create_environment():
    """
    create the tissue environment and initialize the mediator fields.

    returns:
        c: pro-inflammatory mediator concentration matrix.
        g: anti-inflammatory mediator concentration matrix.
    """

    # initialize both chemical concentration fields
    c, g = initialize_fields()

    # initialize the central damaged tissue region
    c = initialize_damage_region(c)

    return c, g

def get_pro_inflammatory(c, x, y):
    """
    return the pro-inflammatory mediator concentration at a given position.
    """

    # apply periodic boundary conditions to the coordinates
    x = x % c.shape[1]
    y = y % c.shape[0]

    # retrieve the local pro-inflammatory concentration
    return c[y, x]


def get_anti_inflammatory(g, x, y):
    """
    return the anti-inflammatory mediator concentration at a given position.
    """

    # apply periodic boundary conditions to the coordinates
    x = x % GRID_SIZE
    y = y % GRID_SIZE

    # retrieve the local anti-inflammatory concentration
    return g[y, x]


def get_concentrations(c, g, x, y):
    """
    return both mediator concentrations at a given position.
    """

    # retrieve the pro-inflammatory concentration
    pro_concentration = get_pro_inflammatory(c, x, y)

    # retrieve the anti-inflammatory concentration
    anti_concentration = get_anti_inflammatory(g, x, y)

    # return both concentrations
    return pro_concentration, anti_concentration


def add_pro_inflammatory(c, x, y, amount):
    """
    increase the pro-inflammatory mediator concentration at a given position.
    """

    # apply periodic boundary conditions
    x = x % GRID_SIZE
    y = y % GRID_SIZE

    # add the specified amount to the local concentration
    c[y, x] += amount


def add_anti_inflammatory(g, x, y, amount):
    """
    increase the anti-inflammatory mediator concentration at a given position.
    """

    # apply periodic boundary conditions
    x = x % GRID_SIZE
    y = y % GRID_SIZE

    # add the specified amount to the local concentration
    g[y, x] += amount


def plot_environment(c, g):
    """
    Display the mediator concentration fields with 100 x 100 gridlines to create the visual environment to validate it
    """

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    fields = [
        (c, "Pro-inflammatory mediator (c)", "Reds"),
        (g, "Anti-inflammatory mediator (g)", "Blues"),
    ]

    for ax, (field, title, cmap) in zip(axes, fields): #we loop through the two matrices for axes, matrix, title, and color scheme of that plot
        #there will be two iterations
        image = ax.imshow(
            field,
            origin="lower",
            cmap=cmap,
            vmin=0,
            vmax=max(field.max(), INITIAL_C),
            interpolation="nearest",
        )
        #this creates an image as 2d array
        
        #label the plot for each iteration
        ax.set_title(title)
        ax.set_xlabel("Grid x-coordinate")
        ax.set_ylabel("Grid y-coordinate")

        #pllace ticks at the boundaries of each grid cell for gridlines
        ax.set_xticks(np.arange(-0.5, GRID_SIZE, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, GRID_SIZE, 1), minor=True)

        # draw gridlines at every cell boundary
        ax.grid(
            which="minor",
            color="black",
            linestyle="-",
            linewidth=0.25,
        )

        # label every 10th coordinate to keep the axes readable
        ax.set_xticks(np.arange(0, GRID_SIZE, 10))
        ax.set_yticks(np.arange(0, GRID_SIZE, 10))

        fig.colorbar(image, ax=ax, label="Concentration")

    plt.tight_layout()
    plt.show()

def calculate_laplacian(field):
    """
    calculate the five-point discrete laplacian of a chemical field.

    the laplacian measures the difference between the concentration
    at a grid location and the concentrations at its four neighbors.
    """

    # retrieve the concentration values immediately to the left
    left = np.roll(field, 1, axis=1)

    # retrieve the concentration values immediately to the right
    right = np.roll(field, -1, axis=1)

    # retrieve the concentration values immediately above
    up = np.roll(field, -1, axis=0)

    # retrieve the concentration values immediately below
    down = np.roll(field, 1, axis=0)

    # calculate the five-point discrete laplacian
    # this measures the local spatial difference in concentration
    laplacian = left + right + up + down - 4 * field

    return laplacian

def update_pro_inflammatory(c, dt, diffusion_c, decay_c):
    """
    update the pro-inflammatory mediator field through diffusion and decay.

    arguments:
        c: current pro-inflammatory concentration matrix.
        dt: numerical timestep.
        diffusion_c: diffusion coefficient for c.
        decay_c: decay rate for c.

    returns:
        updated pro-inflammatory concentration matrix.
    """

    # calculate the laplacian to determine how c diffuses
    laplacian_c = calculate_laplacian(c)

    # calculate the change in concentration caused by diffusion and decay
    change_c = diffusion_c * laplacian_c - decay_c * c

    # update the concentration using the explicit euler method
    c = c + dt * change_c

    # prevent negative concentrations caused by numerical errors
    c = np.maximum(c, 0)

    return c

def update_anti_inflammatory(g, dt, diffusion_g, decay_g):
    """
    update the anti-inflammatory mediator field through diffusion and decay.

    arguments:
        g: current anti-inflammatory concentration matrix.
        dt: numerical timestep.
        diffusion_g: diffusion coefficient for g.
        decay_g: decay rate for g.

    returns:
        updated anti-inflammatory concentration matrix.
    """

    # calculate the laplacian to determine how g diffuses
    laplacian_g = calculate_laplacian(g)

    # calculate the change in concentration caused by diffusion and decay
    change_g = diffusion_g * laplacian_g - decay_g * g

    # update the concentration using the explicit euler method
    g = g + dt * change_g

    # prevent negative concentrations caused by numerical errors
    g = np.maximum(g, 0)

    return g

def update_fields(c, g, dt, diffusion_c, diffusion_g, decay_c, decay_g):
    """
    advance both chemical concentration fields by one timestep.
    """

    # update the pro-inflammatory mediator field
    c = update_pro_inflammatory(c, dt, diffusion_c, decay_c)

    # update the anti-inflammatory mediator field
    g = update_anti_inflammatory(g, dt, diffusion_g, decay_g)

    # return both updated concentration fields
    return c, g

def get_neighborhood_concentration(field, x, y):
    """
    retrieve the chemical concentrations in the eight neighboring
    grid locations surrounding a given position.

    arguments:
        field: chemical concentration matrix (c or g).
        x: current x-coordinate.
        y: current y-coordinate.

    returns:
        a list containing the concentrations at the eight neighboring locations.
    """

    # create an empty list to store neighboring concentrations
    concentrations = []

    # examine the eight surrounding grid locations
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:

            # skip the current location
            if dx == 0 and dy == 0:
                continue

            # calculate the neighboring coordinates
            neighbor_x = (x + dx) % GRID_SIZE
            neighbor_y = (y + dy) % GRID_SIZE

            # retrieve the chemical concentration at this location
            concentration = field[neighbor_y, neighbor_x]

            # add the concentration to the list
            concentrations.append(concentration)

    # return the concentrations of all eight neighbors
    return concentrations

if __name__ == "__main__":

    c, g = create_environment()

    # test chemical concentration retrieval and modification

    # retrieve the initial concentrations at the center of the grid
    pro, anti = get_concentrations(c, g, 50, 50)

    print("\ninitial concentrations at center:")
    print("pro-inflammatory:", pro)
    print("anti-inflammatory:", anti)

    # increase the pro-inflammatory concentration
    add_pro_inflammatory(c, 50, 50, 0.5)

    # increase the anti-inflammatory concentration
    add_anti_inflammatory(g, 50, 50, 0.3)

    # retrieve the updated concentrations
    pro, anti = get_concentrations(c, g, 50, 50)

    print("\nupdated concentrations at center:")
    print("pro-inflammatory:", pro)
    print("anti-inflammatory:", anti)



    # test diffusion and decay

    # initialize the chemical concentration fields
    c, g = create_environment()

    # define temporary numerical parameters for testing
    dt = 0.1
    diffusion_c = 0.1
    diffusion_g = 0.1
    decay_c = 0.01
    decay_g = 0.01

    # record the initial maximum concentration
    initial_max_c = c.max()

    # update both chemical fields for one timestep
    c, g = update_fields(
        c, g, dt,
        diffusion_c, diffusion_g,
        decay_c, decay_g
    )

    # print the results
    print("\nafter one timestep:")
    print("maximum c concentration:", c.max())
    print("maximum g concentration:", g.max())
    print("total c concentration:", np.sum(c))
    print("total g concentration:", np.sum(g))


    # test neighborhood concentration retrieval

    # retrieve the eight neighboring pro-inflammatory concentrations
    neighboring_c = get_neighborhood_concentration(c, 50, 50)

    # print the concentrations
    print("\nneighboring pro-inflammatory concentrations:")
    print(neighboring_c)

    # calculate the average concentration in the neighborhood
    average_c = np.mean(neighboring_c)

    print("average neighboring concentration:", average_c)