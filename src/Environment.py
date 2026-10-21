import numpy as np
import matplotlib.pyplot as plt


# Grid parameters as given in the paper
GRID_SIZE = 100
DAMAGE_RADIUS = 10
INITIAL_C = 1.0


def create_environment():
    """
    Create the tissue grid and initialize the mediator fields.

    Returns:
        c: Pro-inflammatory mediator concentration.
        g: Anti-inflammatory mediator concentration.
    """

    # creating the two chemical concentration fields each as a matrix of zeros
    c = np.zeros((GRID_SIZE, GRID_SIZE), dtype=float)
    g = np.zeros((GRID_SIZE, GRID_SIZE), dtype=float)

    # finding the center of the grid to create initial damage region
    center_x = GRID_SIZE // 2
    center_y = GRID_SIZE // 2

    # creating coordinates for every grid location
    # ogrid returns a multidimensionsal grid using index notation
    x, y = np.ogrid[:GRID_SIZE, :GRID_SIZE] 

    # iddentifying the circular damage region using distance formula
    distance_squared = ((x - center_x) ** 2 + (y - center_y) ** 2)

    damage_region = distance_squared <= DAMAGE_RADIUS ** 2

    # initialize the pro-inflammatory mediator in the damaged region in the array
    c[damage_region] = INITIAL_C 

    return c, g

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
            vmax=INITIAL_C,
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

if __name__ == "__main__":
    c, g = create_environment()

    print("Environment created successfully!")
    print(f"Grid dimensions: {c.shape}")
    print(f"Maximum c concentration: {c.max()}")
    print(f"Maximum g concentration: {g.max()}")
    print(f"Number of damaged grid locations: {np.count_nonzero(c)}")

    plot_environment(c, g)


