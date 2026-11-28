import streamlit as st
import plotly.graph_objects as go

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from simulation import simulation

from visualizations import (
    render_population_metrics,
    render_spatial_visualization,
    render_inflammation_graph,
    render_inflammation_state, 
    render_recruitment_movement_timeline
)

from controls import (
    initialize_simulation,
    render_controls
)


# configure the dashboard
st.set_page_config(
    page_title="Inflammatory Simulation Dashboard",
    page_icon="🧬",
    layout="wide"
)


# define the biological parameters used by the simulation
def create_simulation(pnr=0.02, pmr=0.04):

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
    sim = simulation(
        grid_size=100,
        neutrophil_positions=[],
        macrophage_positions=[],
        neutrophil_params=neutrophil_params,
        macrophage_params=macrophage_params,
        pnr=pnr,
        pmr=pmr
    )

    return sim


# Initialize the simulation
initialize_simulation(create_simulation)


# Dashboard title
st.title("Modeling the Inflammatory Response")

st.write(
    "An interactive spatial simulation of the cellular and molecular "
    "dynamics involved in inflammation."
)


# Create dashboard tabs
background_tab, simulation_tab, discussion_tab = st.tabs(
    [
        "Background & Model",
        "Interactive Simulation",
        "Discussion & Takeaways"
    ]
)


# ------------------------------------------
# TAB 1: BACKGROUND & MODEL
# ------------------------------------------

with background_tab:

    st.header("Understanding Inflammation in Rheumatoid Arthritis")

    st.markdown("""
    Rheumatoid arthritis (RA) is a chronic inflammatory disease in which the
    immune system mistakenly attacks the body's own joint tissues.

    Inflammation is normally a protective response. Immune cells are recruited
    to a site of injury or infection, help remove harmful material, and
    eventually contribute to resolving the response. However, when
    pro-inflammatory signaling persists or the mechanisms responsible for
    resolution are insufficient, inflammation can become chronic.
    """)

    st.info(
        "**Our driving question:** How do changes in neutrophil and macrophage "
        "recruitment influence the balance between pro-inflammatory and "
        "anti-inflammatory signaling over time?"
    )

    st.divider()

    # --------------------------------------------------
    # 1. THE TWO INFLAMMATORY MEDIATORS
    # --------------------------------------------------

    st.subheader("1. Understanding the inflammatory environment")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        <div style="
            background-color:#E8F1FC;
            padding:22px;
            border-radius:12px;
            border-left:5px solid #2563EB;
            min-height:180px;
        ">
            <h3 style="color:#1D4ED8;">
                Pro-inflammatory mediator (c)
            </h3>
            <p style="color:#1F2937; font-size:16px; line-height:1.6;">
                Represents pro-inflammatory signaling in the environment.
                Neutrophils contribute to this signal by releasing
                pro-inflammatory mediator.
            </p>
        </div>
        """, unsafe_allow_html=True)

    with col2:
        st.markdown("""
        <div style="
            background-color:#E8F5EC;
            padding:22px;
            border-radius:12px;
            border-left:5px solid #16A34A;
            min-height:180px;
        ">
            <h3 style="color:#15803D;">
                Anti-inflammatory mediator (g)
            </h3>
            <p style="color:#1F2937; font-size:16px; line-height:1.6;">
                Represents anti-inflammatory signaling. Activated
                macrophages can release this mediator and contribute
                to resolving inflammation.
            </p>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("""
    The simulation tracks both mediators across a spatial grid. Their
    concentrations change over time as cells move, interact, release
    mediators, and are recruited into the environment.
    """)

    st.divider()

    # --------------------------------------------------
    # 2. NEUTROPHIL RECRUITMENT PATHWAY
    # --------------------------------------------------

    st.subheader("2. Neutrophil recruitment pathway")

    st.markdown("""
    Neutrophil recruitment depends on the local inflammatory environment.
    When the recruitment conditions are met, the model uses PNR
    (neutrophil recruitment probability) to determine whether a new
    neutrophil is recruited.
    """)

    neutrophil_pathway = """
    digraph NeutrophilPathway {
        graph [bgcolor="transparent", rankdir=LR, pad="0.25",
               nodesep="0.35", ranksep="0.5"];
        node [shape=box, style="rounded,filled", fontname="Arial",
              fontsize=11, margin="0.16,0.10", color="#64748b"];
        edge [fontname="Arial", fontsize=9, color="#64748b"];

        Environment [label="Local environment\\n(c and g)", fillcolor="#E8F1FC"];
        Conditions [label="Recruitment conditions met?\\nc > 0.05 and g < 0.015",
                    fillcolor="#F1F5F9"];
        PNR [label="PNR\\nRecruitment probability", fillcolor="#DBEAFE"];
        Neutrophil [label="Neutrophil recruited", fillcolor="#DCEFE3"];
        Release [label="Pro-inflammatory\\nmediator release", fillcolor="#FCE8E8"];

        Environment -> Conditions;
        Conditions -> PNR [label="Yes"];
        PNR -> Neutrophil [label="Successful recruitment"];
        Neutrophil -> Release;
        Release -> Environment [label="Changes local c"];
    }
    """

    st.graphviz_chart(neutrophil_pathway, use_container_width=True)

    st.caption(
        "The diagram represents the model's simplified recruitment feedback. "
        "PNR controls recruitment probability; it is not the same as the "
        "probability that a neutrophil releases mediator."
    )

    st.divider()

    # --------------------------------------------------
    # 3. RECRUITMENT AND MOVEMENT SCHEDULE
    # --------------------------------------------------

    st.subheader("3. Recruitment and movement schedule")

    st.markdown("""
    Neutrophils and macrophages operate on different recruitment and
    movement schedules. The timeline below shows when the model checks
    for recruitment and how frequently each cell type moves.
    """)

    render_recruitment_movement_timeline()

    st.divider()

    # --------------------------------------------------
    # 4. MACROPHAGES AND INFLAMMATION RESOLUTION
    # --------------------------------------------------

    st.subheader("4. Macrophages and inflammation resolution")

    st.markdown("""
    Macrophages can contribute to inflammation resolution by engulfing
    apoptotic neutrophils. In this model, a macrophage becomes activated
    after successful phagocytosis. Activated macrophages can then release
    anti-inflammatory mediator (`g`).

    This creates a second pathway that can counterbalance
    pro-inflammatory signaling:
    """)

    macrophage_pathway = """
    digraph MacrophagePathway {
        graph [bgcolor="transparent", rankdir=LR, pad="0.25",
               nodesep="0.35", ranksep="0.5"];
        node [shape=box, style="rounded,filled", fontname="Arial",
              fontsize=11, margin="0.16,0.10", color="#64748b"];
        edge [fontname="Arial", fontsize=9, color="#64748b"];

        Apoptotic [label="Apoptotic\\nneutrophil", fillcolor="#F1F5F9"];
        Phagocytosis [label="Macrophage\\nphagocytosis", fillcolor="#E5F4E9"];
        Activation [label="Macrophage\\nactivation", fillcolor="#DCEFE3"];
        Release [label="Anti-inflammatory\\nmediator release (g)", fillcolor="#E8F1FC"];
        Environment [label="Inflammatory\\nenvironment", fillcolor="#F1F5F9"];

        Apoptotic -> Phagocytosis;
        Phagocytosis -> Activation;
        Activation -> Release;
        Release -> Environment;
    }
    """

    st.graphviz_chart(macrophage_pathway, use_container_width=True)

    st.caption(
        "Macrophages are recruited where the pro-inflammatory mediator "
        "exceeds the model's recruitment threshold (c > 0.4). Their "
        "anti-inflammatory release depends on activation."
    )

    st.divider()

    # --------------------------------------------------
    # 5. FOLLOWING INFLAMMATION OVER TIME
    # --------------------------------------------------

    st.subheader("5. Following inflammation over time")

    st.markdown("""
    The time-series graph in the Interactive Simulation tab tracks the
    mean concentrations of `c` and `g` across the grid at each simulation
    tick. Viewing both signals together helps us investigate how the
    inflammatory environment changes over time.
    """)

    state_col1, state_col2, state_col3 = st.columns(3)

    with state_col1:
        st.markdown("""
        <div style="border-top:5px solid #22C55E; padding:12px;">
            <h4>Stable / Resolved</h4>
            <p>Inflammatory signaling remains relatively low or settles
            toward a stable level.</p>
        </div>
        """, unsafe_allow_html=True)

    with state_col2:
        st.markdown("""
        <div style="border-top:5px solid #F59E0B; padding:12px;">
            <h4>Oscillatory</h4>
            <p>Mediator concentrations fluctuate over time rather than
            settling into a steady pattern.</p>
        </div>
        """, unsafe_allow_html=True)

    with state_col3:
        st.markdown("""
        <div style="border-top:5px solid #EF4444; padding:12px;">
            <h4>Chronic</h4>
            <p>Pro-inflammatory signaling remains elevated, suggesting
            that inflammation is not resolving.</p>
        </div>
        """, unsafe_allow_html=True)

    st.caption(
        "These are the qualitative states we aim to investigate. "
        "The dashboard's current state-classification thresholds are "
        "provisional and should be interpreted as model indicators, "
        "not clinical diagnoses."
    )

    st.divider()


# ------------------------------------------
# TAB 2: INTERACTIVE SIMULATION
# ------------------------------------------

with simulation_tab:

    st.header("Interactive Simulation")

    control_column, visualization_column = st.columns([1, 4])

    with control_column:

        sim = render_controls(create_simulation)

        st.metric("Simulation time", f"{sim.time} ticks")

    with visualization_column:

        # Population metrics
        render_population_metrics(sim)

        st.divider()

        # Spatial map and time-series graph
        spatial_column, graph_column = st.columns(2)

        with spatial_column:

            render_spatial_visualization(sim)

        with graph_column:

            render_inflammation_state(
                st.session_state.inflammation_history
            )

            render_inflammation_graph(
                st.session_state.inflammation_history
            )

            


# ------------------------------------------
# TAB 3: DISCUSSION & TAKEAWAYS
# ------------------------------------------

with discussion_tab:

    st.header("Discussion & Takeaways")

    st.markdown("""
    Explore how immune cell recruitment influences the balance between
    pro-inflammatory and anti-inflammatory signaling. By adjusting
    neutrophil recruitment (PNR) and macrophage recruitment (PMR), we
    can investigate how different conditions shape inflammatory dynamics.
    """)

    st.divider()

    # --------------------------------------------------
    # 1. EXPERIMENTAL COMPARISONS
    # --------------------------------------------------

    st.subheader("1. Experimental Comparisons")

    st.markdown("""
    The simulation allows us to investigate three questions about
    immune cell recruitment and inflammatory behavior.
    """)

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        <div style="
            background:#E8F1FC;
            border-left:5px solid #2563EB;
            border-radius:12px;
            padding:20px;
            min-height:230px;
        ">
            <h3 style="color:#1D4ED8;">01. Neutrophil Recruitment</h3>
            <p style="color:#1F2937;">
                <b>Variable: PNR</b>
            </p>
            <p style="color:#1F2937; line-height:1.6;">
                Compare simulations with different PNR values while
                keeping PMR constant. Examine how changes in neutrophil
                recruitment affect pro-inflammatory mediator levels
                and the overall inflammatory response.
            </p>
        </div>
        """, unsafe_allow_html=True)

    with col2:
        st.markdown("""
        <div style="
            background:#E8F5EC;
            border-left:5px solid #16A34A;
            border-radius:12px;
            padding:20px;
            min-height:230px;
        ">
            <h3 style="color:#15803D;">02. Macrophage Recruitment</h3>
            <p style="color:#1F2937;">
                <b>Variable: PMR</b>
            </p>
            <p style="color:#1F2937; line-height:1.6;">
                Compare simulations with different PMR values while
                keeping PNR constant. Examine how macrophage recruitment
                affects apoptotic neutrophil clearance and
                anti-inflammatory signaling.
            </p>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("<div style='height:12px'></div>", unsafe_allow_html=True)

    st.markdown("""
    <div style="
        background:#F0F5FA;
        border-left:5px solid #64748B;
        border-radius:12px;
        padding:20px;
    ">
        <h3 style="color:#334155;">03. Combined Recruitment Effects</h3>
        <p style="color:#1F2937;">
            <b>Variables: PNR and PMR</b>
        </p>
        <p style="color:#1F2937; line-height:1.6;">
            Explore how different combinations of neutrophil and
            macrophage recruitment influence the interaction between
            pro-inflammatory and anti-inflammatory mediators.
        </p>
    </div>
    """, unsafe_allow_html=True)

    st.divider()

    # --------------------------------------------------
    # 2. INFLAMMATORY RESPONSE PATTERNS
    # --------------------------------------------------

    st.subheader("2. Inflammatory Response Patterns")

    st.markdown("""
    Changes in mediator concentrations over time can be used to
    investigate three possible patterns of inflammatory behavior.
    """)

    state1, state2, state3 = st.columns(3)

    with state1:
        st.markdown("""
        <div style="
            background:#E8F5EC;
            border-top:5px solid #16A34A;
            border-radius:12px;
            padding:18px;
            min-height:210px;
        ">
            <h3 style="color:#15803D;">Stable / Resolved</h3>
            <p style="color:#1F2937; line-height:1.6;">
                Inflammatory activity remains low or decreases over
                time, suggesting that the response is being regulated.
            </p>
        </div>
        """, unsafe_allow_html=True)

    with state2:
        st.markdown("""
        <div style="
            background:#E8F1FC;
            border-top:5px solid #2563EB;
            border-radius:12px;
            padding:18px;
            min-height:210px;
        ">
            <h3 style="color:#1D4ED8;">Oscillatory</h3>
            <p style="color:#1F2937; line-height:1.6;">
                Inflammatory mediator levels fluctuate over time,
                suggesting repeated changes in inflammatory activity.
            </p>
        </div>
        """, unsafe_allow_html=True)

    with state3:
        st.markdown("""
        <div style="
            background:#F0F5FA;
            border-top:5px solid #475569;
            border-radius:12px;
            padding:18px;
            min-height:210px;
        ">
            <h3 style="color:#334155;">Chronic</h3>
            <p style="color:#1F2937; line-height:1.6;">
                Inflammatory activity remains elevated over time,
                suggesting that the response is not resolving.
            </p>
        </div>
        """, unsafe_allow_html=True)

    st.caption(
        "These are possible patterns to investigate through the simulation, "
        "not predetermined experimental outcomes."
    )

    st.divider()

    # --------------------------------------------------
    # 3. KEY TAKEAWAYS
    # --------------------------------------------------

    st.subheader("3. Key Takeaways")

    st.markdown("""
    <div style="
        background:#E8F1FC;
        border-radius:12px;
        padding:22px;
        border-left:5px solid #2563EB;
    ">
        <h3 style="color:#1D4ED8;">The central question</h3>
        <p style="color:#1F2937; font-size:16px; line-height:1.7;">
            How does the recruitment of different immune cell populations
            influence inflammatory dynamics, and are changes in PNR or
            PMR associated with stable, oscillatory, or chronic
            inflammatory patterns?
        </p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    Through these comparisons, we aim to understand how the interaction
    between neutrophils, macrophages, and inflammatory mediators
    contributes to the regulation or persistence of inflammation.
    """)

    st.divider()

    # --------------------------------------------------
    # 4. LIMITATIONS AND FUTURE DIRECTIONS
    # --------------------------------------------------

    st.subheader("4. Limitations & Future Directions")

    limitation_col, future_col = st.columns(2)

    with limitation_col:
        st.markdown("""
        <div style="
            background:#F1F5F9;
            border-radius:12px;
            padding:20px;
            min-height:230px;
        ">
            <h3 style="color:#334155;">Limitations</h3>
            <ul style="color:#1F2937; line-height:1.8;">
                <li>Simplified representation of inflammation</li>
                <li>Does not capture every immune cell type</li>
                <li>Limited representation of biological signaling pathways</li>
                <li>Model behavior depends on selected parameters</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)

    with future_col:
        st.markdown("""
        <div style="
            background:#E8F5EC;
            border-radius:12px;
            padding:20px;
            min-height:230px;
        ">
            <h3 style="color:#15803D;">Future Directions</h3>
            <ul style="color:#1F2937; line-height:1.8;">
                <li>Explore additional immune cell interactions</li>
                <li>Refine model parameters</li>
                <li>Compare repeated simulation runs</li>
                <li>Investigate consistency across recruitment conditions</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)

    st.divider()

    st.caption(
        "This computational model is a simplified representation of "
        "inflammatory dynamics and is not intended to reproduce every "
        "biological process involved in rheumatoid arthritis."
    )