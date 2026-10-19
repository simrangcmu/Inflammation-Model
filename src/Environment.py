import numpy as np
import matplotlib.pyplot as plt


# Grid parameters
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

    # Create the two chemical concentration fields.
    c = np.zeros((GRID_SIZE, GRID_SIZE), dtype=float)
    g = np.zeros((GRID_SIZE, GRID_SIZE), dtype=float)

    # Find the center of the grid.
    center_x = GRID_SIZE // 2
    center_y = GRID_SIZE // 2

    # Create coordinates for every grid location.
    x, y = np.ogrid[:GRID_SIZE, :GRID_SIZE]

    # Identify the circular damage region.
    distance_squared = (
        (x - center_x) ** 2
        + (y - center_y) ** 2
    )

    damage_region = distance_squared <= DAMAGE_RADIUS ** 2

    # Initialize the pro-inflammatory mediator in the damaged region.
    c[damage_region] = INITIAL_C

    return c, g


def plot_environment(c, g):
    """
    Display the initial pro-inflammatory and anti-inflammatory
    mediator concentration fields.
    """

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Plot the pro-inflammatory mediator.
    image_c = axes[0].imshow(
        c,
        origin="lower",
        cmap="Reds",
        vmin=0,
        vmax=INITIAL_C,
    )
    axes[0].set_title("Pro-inflammatory mediator (c)")
    axes[0].set_xlabel("Grid x-coordinate")
    axes[0].set_ylabel("Grid y-coordinate")
    fig.colorbar(image_c, ax=axes[0], label="Concentration")

    # Plot the anti-inflammatory mediator.
    image_g = axes[1].imshow(
        g,
        origin="lower",
        cmap="Blues",
        vmin=0,
        vmax=INITIAL_C,
    )
    axes[1].set_title("Anti-inflammatory mediator (g)")
    axes[1].set_xlabel("Grid x-coordinate")
    axes[1].set_ylabel("Grid y-coordinate")
    fig.colorbar(image_g, ax=axes[1], label="Concentration")

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
