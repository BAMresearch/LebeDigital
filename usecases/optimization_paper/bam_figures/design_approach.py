import os

# from pathlib import Path
import pathlib

import graphviz
import matplotlib.pyplot as plt

# looks better but different to other text
plt.rc("mathtext", fontset="cm")  # matplotlib: force computer modern font set

LOCATION = pathlib.Path(__file__).parent.resolve()
IMAGE_FOLDER = LOCATION / "graph_tex_nodes"
pathlib.Path(IMAGE_FOLDER).mkdir(exist_ok=True)

FONTSIZE = 12


def design_approach_graph(file_name_1, file_name_2, view=False):
    file_name_1 = str(file_name_1)  # convert pathlib object to useful string
    file_name_2 = str(file_name_2)  # convert pathlib object to useful string

    process = "blue"
    input = "green"
    problem = "red"
    kpi = "orange"

    def create_nodes(file_name="test_graph", header=None):
        dot = graphviz.Digraph(file_name, comment="Design Approach", format="pdf")
        dot.attr(rankdir="LR")
        dot.attr("node", fontsize=str(FONTSIZE))
        if header:
            dot.attr(label=header, labelloc="t", labeljust="c", fontsize=str(FONTSIZE * 1.5))

        dot.node("loads", "Loads")
        dot.node("structural_constraint", "Structural Constraint")
        dot.node("structural_design", "Structural Design", shape="rectangle", color=process)
        dot.node("height", "Beam Height")
        dot.node("mechanical_properties", "Mechanical Properties")
        dot.node("mix_design", "Concrete Mix Design", shape="rectangle", color=process)
        dot.node("mix_constraint", "Mix Constraint")
        dot.node("slag_content", "Slag Content")
        return dot

    dot = create_nodes(file_name_1, "Standard Design Approach")

    dot.edge("loads", "structural_design")
    dot.edge("structural_constraint", "structural_design")
    dot.edge("structural_design", "height")
    dot.edge("structural_design", "mechanical_properties")
    dot.edge("mechanical_properties", "mix_design")
    dot.edge("mix_constraint", "mix_design")
    dot.edge("mix_design", "slag_content")

    dot.render(view=view, cleanup=True)

    dot = create_nodes(file_name_2, "Proposed Workflow")

    dot.edge("loads", "structural_design")
    dot.edge("structural_constraint", "structural_design")
    dot.edge("height", "structural_design")
    dot.edge("mechanical_properties", "structural_design")
    dot.edge("mix_design", "mechanical_properties")
    dot.edge("mix_constraint", "mix_design")
    dot.edge("slag_content", "mix_design")

    dot.render(view=view, cleanup=True)


if __name__ == "__main__":
    design_approach_graph(file_name_1="test_standard_design", file_name_2="test_proposed_design", view=True)
