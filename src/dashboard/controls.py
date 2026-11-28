
import streamlit as st
import numpy as np


def initialize_simulation(create_simulation):
    """Create the initial simulation and initialize its history."""

    if "simulation" not in st.session_state:
        st.session_state.simulation = create_simulation(
            pnr=0.02,
            pmr=0.04
        )

    if "inflammation_history" not in st.session_state:
        reset_history(st.session_state.simulation)


def reset_history(sim):
    """Start a new time series for a simulation."""

    mean_c = float(np.mean(sim.environment[0]))
    mean_g = float(np.mean(sim.environment[1]))

    st.session_state.inflammation_history = {
        "time": [sim.time],
        "mean_c": [mean_c],
        "mean_g": [mean_g]
    }


def record_inflammation(sim):
    """Record average mediator concentrations after a simulation step."""

    history = st.session_state.inflammation_history

    history["time"].append(sim.time)
    history["mean_c"].append(float(np.mean(sim.environment[0])))
    history["mean_g"].append(float(np.mean(sim.environment[1])))


def render_controls(create_simulation):
    """Display simulation controls and return the current simulation."""

    st.subheader("Simulation Parameters")

    pnr = st.slider(
        "PNR — Neutrophil Recruitment Probability",
        min_value=0.0,
        max_value=1.0,
        value=0.02,
        step=0.01,
        help="Probability of recruiting neutrophils at eligible grid locations."
    )

    pmr = st.slider(
        "PMR — Macrophage Recruitment Probability",
        min_value=0.0,
        max_value=1.0,
        value=0.04,
        step=0.01,
        help="Probability of recruiting macrophages at eligible grid locations."
    )

    st.caption(
        "PNR controls neutrophil recruitment. PMR controls macrophage "
        "recruitment. The initial values are based on the paper."
    )

    # Reinitialize if either parameter changes.
    sim = st.session_state.simulation

    if pnr != sim.pnr or pmr != sim.pmr:
        st.session_state.simulation = create_simulation(
            pnr=pnr,
            pmr=pmr
        )
        reset_history(st.session_state.simulation)
        sim = st.session_state.simulation

    st.divider()
    st.subheader("Simulation Controls")

    col1, col2 = st.columns(2)

    with col1:
        if st.button("Advance 1 tick", use_container_width=True):
            sim.step()
            record_inflammation(sim)

    with col2:
        if st.button("Advance 10 ticks", use_container_width=True):
            for _ in range(10):
                sim.step()
                record_inflammation(sim)

    if st.button("Reset simulation", use_container_width=True):
        st.session_state.simulation = create_simulation(
            pnr=pnr,
            pmr=pmr
        )
        reset_history(st.session_state.simulation)
        st.rerun()

    return st.session_state.simulation
