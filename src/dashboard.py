import streamlit as st
import matplotlib.pyplot as plt

from Environment import create_environment, GRID_SIZE, DAMAGE_RADIUS


states = {
    "Resolution": create_environment(),
    "Persistent": create_environment(),
    "Oscillatory": create_environment()
}

#initial state (with initial conditions)
col1, col2, col3 = st.columns(3)

#col1.metric("Grid size", f"{GRID_SIZE} × {GRID_SIZE}")
col2.metric("Damage radius", DAMAGE_RADIUS)
col3.metric(
    "Initial damaged locations",
    int((states["Resolution"][0] > 0).sum())
)

st.divider()

#3 spatial grids (one for each state)

grid_columns = st.columns(3)

for column, (scenario_name, fields) in zip(
    grid_columns, states.items()
):
    c, g = fields

    with column:
        st.markdown(f"### {scenario_name}")

        fig, ax = plt.subplots(figsize=(4, 4))

        image = ax.imshow(
            c,
            origin="lower",
            cmap="Reds",
            vmin=0,
            vmax=1
        )

        ax.set_title("Pro-inflammatory field (c)")
        ax.set_xlabel("X-coordinate")
        ax.set_ylabel("Y-coordinate")

        fig.colorbar(image, ax=ax, label="Concentration")

        st.pyplot(fig)
        plt.close(fig)

        st.caption("Initial state — time step 0")

st.divider()

#time graphs (one for each state) 
graph_columns = st.columns(3)

for column, scenario_name in zip(
    graph_columns, states.keys()
):
    with column:
        st.markdown(f"### {scenario_name}")

        fig, ax = plt.subplots(figsize=(4, 3))

        ax.set_title("Inflammation vs. Time")
        ax.set_xlabel("Time step")
        ax.set_ylabel("Inflammation level")

        ax.set_xlim(0, 100)
        ax.set_ylim(0, 1)

        ax.text(
            0.5, 0.5,
            "Waiting for simulation data...",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes
        )

        st.pyplot(fig)
        plt.close(fig)
