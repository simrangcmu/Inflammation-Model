import numpy as np

def get_neighbors(position, grid_shape):
    """
    Return the eight neighboring locations surrounding a cell

    Parameters:
        position: Current (row, column) position
        grid_shape: Dimensions of the tissue grid

    Returns:
        List of valid neighboring positions
    """

    row, col = position #representign the current position of the cell in the grid
    rows, cols = grid_shape

    neighbors = [] #empty list to store the valid neighboring positions

    # First, we chheck all eight surrounding directions of the position
    #we examine cells immediately surrounding current position
    for dr in [-1, 0, 1]: #for change in row direction
        for dc in [-1, 0, 1]: #change in column direction
            #hence, we have 9 combinations of these changes
            # excluding the current position as we don't need to consider it
            if dr == 0 and dc == 0: #we dont want to include same spot as movement destination
                continue

            new_row = row + dr #we calulate the new row and column indices for the neighboring position
            new_col = col + dc #we find new posiiton

            # Check that the new position is inside the grid or not
            if 0 <= new_row < rows and 0 <= new_col < cols:
                neighbors.append((new_row, new_col)) #add it to list if yes

    return neighbors #all valid neighboring positions

def calculate_movement_probabilities(
    c,
    position,
    previous_angle=None,
    k_grad=80,
    sigma_mem=1.0
):
    """
    Calculate movement probabilities using chemical gradients
    and directional persistence

    Parameters:
        c: Pro-inflammatory mediator concentration matrix
        position: Current (row, column) position of neutrophil
        previous_angle: Direction of the previous movement in radians
        k_grad: Chemotactic sensitivity
        sigma_mem: Directional persistence parameter (memory)

    Returns:
        neighbors: Valid neighboring positions
        probabilities: Probability of moving to each neighbor
    """

    row, col = position

    # Get the eight neighboring locations
    neighbors = get_neighbors(position, c.shape) #passing the current position of neutrophil and the dimensions of the concentration matrix

    # Current mediator concentration
    current_c = c[row, col]

    # Calculate gradient-based weights in empty list
    gradient_weights = []

    for neighbor in neighbors: #looping over all neighbors

        neighbor_c = c[neighbor] #get the concentration of the pro-inflammatory mediator at the neighboring position

        difference = neighbor_c - current_c #conc difference

        weight = np.exp(k_grad * difference) #weight is calculated using the exponential function, which amplifies the effect of the concentration difference based on the chemotactic sensitivity parameter k_grad

        gradient_weights.append(weight) #adds the calculated weight to list

    gradient_weights = np.array(gradient_weights)

    # Normalize gradient weights to convert to probabilities
    gradient_probabilities = (
        gradient_weights / np.sum(gradient_weights)
    )

    # Calculate directional persistence weights by calculating how much memory neutrophil has of previous direction
    memory_weights = []

    for neighbor in neighbors: #store one memory weight for each neighboring cell

        dr = neighbor[0] - row
        dc = neighbor[1] - col

        # Calculate angle of the potential movement
        angle = np.arctan2(dr, dc) #takes two values and calculates an angle in radians

        if previous_angle is None:

            # No previous direction for before: no memory preference, so each direction gets equal memory weight
            memory_weight = 1.0

        else:

            # Calculate the smallest circular angular difference (how closely each potential movement dir matches previosu dir)
            angle_difference = np.arctan2(
                np.sin(angle - previous_angle),
                np.cos(angle - previous_angle)
            )
            #circular because it needs to wrap around to find accurate angle

            #we use a gaussian shaped direction preference to calc memory weight
            #the smaller the angular diff, larger the memory weight
            memory_weight = np.exp(
                -(angle_difference ** 2) / (2 * sigma_mem ** 2)
            )

        memory_weights.append(memory_weight) #add to list

    memory_weights = np.array(memory_weights)

    # Combine gradient and memory probabilities using element-wise multiplication
    combined_weights = gradient_probabilities * memory_weights

    #neighboring cell has high combined weight when its inflamatory conc is attarctive, 
    # and movement dir is consistent with prrevious dir
    # Normalize the final movement probabilities

    probabilities = combined_weights / np.sum(combined_weights)

    return neighbors, probabilities #each prob corresponds to that neighbor

def choose_next_position(
    c,
    position,
    previous_angle=None,
    k_grad=80,
    sigma_mem=1.0
):
    """
    Select the cell's next position using chemotaxis probabilities

    Returns:
        new_position: Selected neighboring location
        new_angle: Direction of the movement in radians
    """

    neighbors, probabilities = calculate_movement_probabilities(
        c,
        position,
        previous_angle,
        k_grad,
        sigma_mem
    )

    # Randomly select a neighbor according to its probability
    selected_index = np.random.choice(
        len(neighbors),
        p=probabilities
    )

    new_position = neighbors[selected_index] #identify selected position for new position of neutrophil

    # Calculate the direction of the selected movement
    dr = new_position[0] - position[0]
    dc = new_position[1] - position[1]

    new_angle = np.arctan2(dr, dc)

    return new_position, new_angle

if __name__ == "__main__":

    # Create a test concentration field for checking
    c = np.zeros((100, 100))

    #create a concentration gradient increasing toward the right
    for col in range(100): #iterates over every column
        c[:, col] = col / 100 #selects every row in the column and gives it a conc value by dividing by 100

    # Place a test neuttophil cell at the center
    position = (50, 50)

    #chemical gradient will favor movmeent toward the right

    # Calculate movement probabilities
    neighbors, probabilities = calculate_movement_probabilities(
        c,
        position,
        previous_angle=None,
        k_grad=80,
        sigma_mem=1.0
    )

    print("Current position:", position)

    print("\nMovement probabilities:")

    for neighbor, probability in zip(neighbors, probabilities):
        print(f"{neighbor}: {probability:.4f}")

    print("\nSum of probabilities:", probabilities.sum()) #adding all the values in the array to check if 1

    # Simulate one movement
    new_position, new_angle = choose_next_position(
        c,
        position,
        previous_angle=None
    )

    print("\nNew position:", new_position)
    print("New angle:", new_angle)