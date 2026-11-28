
from panel import state
import streamlit as st
import plotly.graph_objects as go
import numpy as np

from active_neutrophils import Neutrophil
from apoptotic_neutrophils import ApoptoticNeutrophil
from macrophage import Macrophage


# --------------------------------------------------
# POPULATION METRICS
# --------------------------------------------------

def render_population_metrics(sim):

    active = sim.grid.get_agents_by_type(Neutrophil)

    apoptotic = sim.grid.get_agents_by_type(
        ApoptoticNeutrophil
    )

    macrophages = sim.grid.get_agents_by_type(Macrophage)

    col1, col2, col3 = st.columns(3)

    col1.metric("Active neutrophils", len(active))
    col2.metric("Apoptotic neutrophils", len(apoptotic))
    col3.metric("Macrophages", len(macrophages))


# --------------------------------------------------
# SPATIAL VISUALIZATION
# --------------------------------------------------

def create_spatial_figure(sim, view):

    c = sim.environment[0]
    g = sim.environment[1]

    figure = go.Figure()

    if view == "Overall Inflammation Heatmap":

        net_inflammation = c - g

        figure.add_trace(
            go.Heatmap(
                z=net_inflammation,
                colorscale="RdBu_r",
                colorbar=dict(title="Net inflammation")
            )
        )

    elif view == "Pro-inflammatory Mediator (c)":

        figure.add_trace(
            go.Heatmap(
                z=c,
                colorscale="Reds",
                colorbar=dict(title="Mediator c")
            )
        )

    elif view == "Anti-inflammatory Mediator (g)":

        figure.add_trace(
            go.Heatmap(
                z=g,
                colorscale="Blues",
                colorbar=dict(title="Mediator g")
            )
        )

    figure.update_layout(
        title=view,
        xaxis_title="X position",
        yaxis_title="Y position",
        height=650,
        margin={"l": 20, "r": 20, "t": 55, "b": 50}
    )

    figure.update_yaxes(
        scaleanchor="x",
        scaleratio=1
    )

    return figure


def render_spatial_visualization(sim):

    st.subheader("Spatial Distribution")

    view = st.selectbox(
        "Select visualization",
        [
            "Overall Inflammation Heatmap",
            "Pro-inflammatory Mediator (c)",
            "Anti-inflammatory Mediator (g)"
        ]
    )

    figure = create_spatial_figure(sim, view)

    st.plotly_chart(
        figure,
        use_container_width=True
    )


# --------------------------------------------------
# INFLAMMATION TIME SERIES
# --------------------------------------------------

def create_inflammation_figure(history):

    figure = go.Figure()

    times = history["time"]
    mean_c = history["mean_c"]
    mean_g = history["mean_g"]

    # Pro-inflammatory mediator
    figure.add_trace(
        go.Scatter(
            x=times,
            y=mean_c,
            mode="lines",
            name="Pro-inflammatory c",
            line={"color": "firebrick", "width": 2.5},
            hovertemplate=(
                "Time: %{x}<br>"
                "Mean c: %{y:.5f}"
                "<extra></extra>"
            )
        )
    )

    # Anti-inflammatory mediator
    figure.add_trace(
        go.Scatter(
            x=times,
            y=mean_g,
            mode="lines",
            name="Anti-inflammatory g",
            line={"color": "royalblue", "width": 2.5},
            hovertemplate=(
                "Time: %{x}<br>"
                "Mean g: %{y:.5f}"
                "<extra></extra>"
            )
        )
    )

    # Add hoverable lifespan labels directly on the graph.
    # These are informational markers, not simulation data.
    if times:

        x_start = times[0]
        x_end = times[-1]

        x_range = max(x_end - x_start, 1)

        current_max = max(
            max(mean_c, default=0),
            max(mean_g, default=0),
            0.001
        )

        label_y = current_max * 1.18

        figure.add_trace(
            go.Scatter(
                x=[
                    x_start + x_range * 0.20,
                    x_start + x_range * 0.65
                ],
                y=[label_y, label_y],
                mode="markers+text",
                text=[
                    "Neutrophil lifespan ⓘ",
                    "Macrophage lifespan ⓘ"
                ],
                textposition="top center",
                marker={
                    "size": 9,
                    "color": ["firebrick", "royalblue"],
                    "symbol": "circle"
                },
                name="Cell lifespan information",
                showlegend=False,
                customdata=[
                    ["60–1,440 simulation ticks"],
                    ["24 hours–60 days (paper reference range)"]
                ],
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "%{customdata[0]}"
                    "<extra></extra>"
                )
            )
        )

        figure.update_yaxes(
            range=[0, label_y * 1.25]
        )

    figure.update_layout(
        title="Inflammatory Mediator Concentrations Over Time",
        xaxis_title="Simulation time (ticks)",
        yaxis_title="Mean mediator concentration",
        height=450,
        hovermode="closest",
        legend={"orientation": "h"},
        margin={"l": 50, "r": 30, "t": 75, "b": 50}
    )

    figure.update_xaxes(
        tickmode="linear",
        tick0=0,
        dtick=1
    )

        # Improve spacing and prevent overlapping labels
    figure.update_layout(
        margin=dict(l=70, r=30, t=100, b=100),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.12,
            xanchor="center",
            x=0.5
        )
    )

    # Show fewer x-axis labels while keeping tick values as integers
    figure.update_xaxes(
        title_text="Simulation time (ticks)",
        tickmode="auto",
        nticks=8,
        tickformat="d",
        automargin=True
    )

    return figure


def render_inflammation_graph(history):

    st.subheader("Inflammatory Dynamics")

    figure = create_inflammation_figure(history)

    st.plotly_chart(
        figure,
        use_container_width=True
    )

    st.caption(
        "Hover over the lifespan labels on the graph for cell lifespan "
        "information. Concentration curves show the mean c and g values "
        "across the environment."
    )


# --------------------------------------------------
# INFLAMMATION STATE
# --------------------------------------------------

def render_inflammation_state(history):
    st.subheader("Inflammation State")

    c_values = history["mean_c"]
    g_values = history["mean_g"]

    if len(c_values) < 20:
        st.info("Run the simulation for more ticks to determine the inflammation state.")
        return

    # Examine the most recent 20 simulation observations
    recent_c = np.array(c_values[-20:])
    recent_g = np.array(g_values[-20:])

    # Calculate the average concentrations
    mean_c = np.mean(recent_c)
    mean_g = np.mean(recent_g)

    # Measure how much pro-inflammatory concentration changes
    variation_c = np.std(recent_c)

    # Determine the current state
    if variation_c > 0.1 * abs(mean_c):
        state = "Oscillatory"
    elif mean_c > 0.4 and mean_g < mean_c:
        state = "Chronic"
    else:
        state = "Stable / Resolved"

    # Assign a color to each state
    state_colors = {
        "Stable / Resolved": "#22C55E",
        "Oscillatory": "#F59E0B",
        "Chronic": "#EF4444"
    }

    color = state_colors[state]

    st.markdown(
        f"""
        <div style="display:flex; align-items:center; gap:9px;">
            <span style="
                height:12px;
                width:12px;
                background-color:{color};
                border-radius:50%;
                display:inline-block;
            "></span>
            <span style="font-size:16px; font-weight:600;">
                {state}
            </span>
        </div>
        """,
        unsafe_allow_html=True
    )

def render_recruitment_movement_timeline():

    st.markdown("**Recruitment and movement schedule**")

    # Simulation ticks shown on the timeline
    ticks = list(range(14))

    # Recruitment occurs at these ticks
    neutrophil_recruitment = list(range(0, 14, 2))
    macrophage_recruitment = list(range(0, 14, 5))

    # Smooth curves illustrating movement frequency
    x = np.linspace(0, 13, 1000)

    # Neutrophils move twice per tick
    neutrophil_movement = (
        1.15 + 0.18 * np.sin(2 * np.pi * 2 * x)
    )

    # Macrophages move once per tick
    macrophage_movement = (
        0.45 + 0.18 * np.sin(2 * np.pi * x)
    )

    figure = go.Figure()

    # Neutrophil recruitment markers
    figure.add_trace(go.Scatter(
        x=neutrophil_recruitment,
        y=[2.65] * len(neutrophil_recruitment),
        mode="markers",
        name="Neutrophil recruitment",
        marker=dict(
            color="#4385C5",
            size=12,
            symbol="circle"
        )
    ))

    # Macrophage recruitment markers
    figure.add_trace(go.Scatter(
        x=macrophage_recruitment,
        y=[2.15] * len(macrophage_recruitment),
        mode="markers",
        name="Macrophage recruitment",
        marker=dict(
            color="#22A06B",
            size=12,
            symbol="circle"
        )
    ))

    # Neutrophil movement curve
    figure.add_trace(go.Scatter(
        x=x,
        y=neutrophil_movement,
        mode="lines",
        name="Neutrophil movement: 2 times/tick",
        line=dict(color="#4385C5", width=2)
    ))

    # Macrophage movement curve
    figure.add_trace(go.Scatter(
        x=x,
        y=macrophage_movement,
        mode="lines",
        name="Macrophage movement: 1 time/tick",
        line=dict(color="#22A06B", width=2)
    ))

    figure.update_layout(
        title="Recruitment Events and Movement Frequency",
        height=480,
        margin=dict(l=175, r=25, t=65, b=65),
        xaxis=dict(
            title="Simulation time (ticks)",
            tickmode="array",
            tickvals=ticks,
            ticktext=[str(tick) for tick in ticks],
            range=[-0.3, 13.3],
            showgrid=True,
            gridcolor="rgba(128,128,128,0.2)"
        ),
        yaxis=dict(
            tickmode="array",
            tickvals=[2.65, 2.15, 1.15, 0.45],
            ticktext=[
                "Neutrophil recruitment",
                "Macrophage recruitment",
                "Neutrophil movement",
                "Macrophage movement"
            ],
            range=[0, 3],
            showgrid=False,
            zeroline=False
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5
        )
    )

    st.plotly_chart(figure, use_container_width=True)

    st.caption(
        "Recruitment markers show the ticks at which the model checks "
        "for new cells. The curves illustrate movement frequency: "
        "neutrophils move twice per tick, while macrophages move once. "
        "The curves are conceptual representations, not measured cell paths."
    )

