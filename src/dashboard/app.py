import streamlit as st
import plotly.graph_objects as go

from simulation import Simulation
from active_neutrophils import Neutrophil
from apoptotic_neutrophils import ApoptoticNeutrophil
from macrophage import Macrophage


# configure the dashboard
st.set_page_config(
    page_title="Inflammatory Simulation Dashboard",
    page_icon="🧬",
    layout="wide"
)


# define the biological parameters used by the simulation
def create_simulation():

    # neutrophil parameters
    neutrophil_params = {
        "movesPerTick_n": 2,
        "n_k_grad": 80,
        "n_sigma": 1.0,
        "prob_neutroReleasePro": 0.5,
        "inc_neutroReleasePro": 0.001
    }

    # macrophage parameters
    macrophage_params = {
        "movesPerTick_m": 1,
        "prob_macroPhagoApop": 1.0,
        "prob_macroReleaseAnti": 0.8,
        "inc_macroReleaseAnti": 0.001,
        "m_k_grad": 80,
        "m_sigma": 1.0,
        "prob_macroLeave": 0.8,
        "thresh_macroLeave": 0.02
    }

    # create a new simulation with no initial agents
    simulation = Simulation(
        grid_size=100,
        neutrophil_positions=[],
        macrophage_positions=[],
        neutrophil_params=neutrophil_params,
        macrophage_params=macrophage_params
    )

    return simulation


# initialize the simulation when the dashboard first opens
if "simulation" not in st.session_state:
    st.session_state.simulation = create_simulation()


# retrieve the current simulation
simulation = st.session_state.simulation


# dashboard title
st.title("Modeling the Inflammatory Response")

st.write(
    "An interactive spatial simulation of the cellular and molecular "
    "dynamics involved in inflammation."
)


# simulation controls
st.subheader("Simulation Controls")

col1, col2, col3 = st.columns(3)

with col1:
    if st.button("Advance 1 tick"):
        simulation.step()

with col2:
    if st.button("Advance 10 ticks"):
        for step_number in range(10):
            simulation.step()

with col3:
    if st.button("Reset simulation"):
        st.session_state.simulation = create_simulation()
        st.rerun()


# display simulation time
st.metric(
    "Simulation time",
    f"{simulation.time} minutes"
)


# retrieve the current populations
active_neutrophils = simulation.grid.get_agents_by_type(Neutrophil)

apoptotic_neutrophils = simulation.grid.get_agents_by_type(
    ApoptoticNeutrophil
)

macrophages = simulation.grid.get_agents_by_type(Macrophage)


# display population counts
st.subheader("Cell Populations")

col1, col2, col3 = st.columns(3)

col1.metric("Active neutrophils", len(active_neutrophils))

col2.metric("Apoptotic neutrophils", len(apoptotic_neutrophils))

col3.metric("Macrophages", len(macrophages))


# retrieve the current mediator fields
c = simulation.environment[0]
g = simulation.environment[1]


# create a function to display the spatial environment and cells
def create_spatial_figure(field, title, colorscale, simulation):

    figure = go.Figure()

    # display the mediator concentration as a heatmap
    figure.add_trace(
        go.Heatmap(
            z=field,
            colorscale=colorscale,
            colorbar={"title": "Concentration"},
            name=title
        )
    )

    # get the current positions of each cell type
    active_positions = [
        simulation.grid.get_location(agent)
        for agent in simulation.grid.get_agents_by_type(Neutrophil)
    ]

    apoptotic_positions = [
        simulation.grid.get_location(agent)
        for agent in simulation.grid.get_agents_by_type(
            ApoptoticNeutrophil
        )
    ]

    macrophage_positions = [
        simulation.grid.get_location(agent)
        for agent in simulation.grid.get_agents_by_type(Macrophage)
    ]

    # display active neutrophils
    if active_positions:
        figure.add_trace(
            go.Scatter(
                x=[position[0] for position in active_positions],
                y=[position[1] for position in active_positions],
                mode="markers",
                marker={"color": "lime", "size": 6, "symbol": "circle"},
                name="Active neutrophils"
            )
        )

    # display apoptotic neutrophils
    if apoptotic_positions:
        figure.add_trace(
            go.Scatter(
                x=[position[0] for position in apoptotic_positions],
                y=[position[1] for position in apoptotic_positions],
                mode="markers",
                marker={"color": "yellow", "size": 6, "symbol": "x"},
                name="Apoptotic neutrophils"
            )
        )

    # display macrophages
    if macrophage_positions:
        figure.add_trace(
            go.Scatter(
                x=[position[0] for position in macrophage_positions],
                y=[position[1] for position in macrophage_positions],
                mode="markers",
                marker={"color": "cyan", "size": 8, "symbol": "diamond"},
                name="Macrophages"
            )
        )

    # configure the spatial display
    figure.update_layout(
        title=title,
        xaxis_title="X position",
        yaxis_title="Y position",
        height=550,
        legend={"orientation": "h"}
    )

    figure.update_yaxes(
        scaleanchor="x",
        scaleratio=1
    )

    return figure


# display both mediator fields
st.subheader("Spatial Distribution of Mediators and Cells")

left_column, right_column = st.columns(2)

with left_column:

    figure_c = create_spatial_figure(
        c,
        "Pro-inflammatory mediator (c)",
        "Reds",
        simulation
    )

    st.plotly_chart(figure_c, use_container_width=True)

with right_column:

    figure_g = create_spatial_figure(
        g,
        "Anti-inflammatory mediator (g)",
        "Blues",
        simulation
    )

    st.plotly_chart(figure_g, use_container_width=True)