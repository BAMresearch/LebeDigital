## Overview
This directory contains to code to automatically generate the planned paper.
This includes the scripts to generate the figures and at least the option to re-run the computations if needed.

## Structure
The current proposal is to write the paper in the **tex** folder and split the document in sections, which are located in **tex/tex_sections**.
Modules like **workflow_graph** can be added as directory to the main directory and the included in the pydoit workflow.
Figures generated for the paper are located in **figures**.

## Features

### pydoit
There is a **dodo.py** file which includes all steps required to generate the paper.
It has the option for different modes to be implemented.
As example the default is currently
```
doit mode=default
```
If there are expensive computations that should be skipped, or settings changed for integration testing purposes, this can be done using this option.

### Easy list of figure dependencies
A function has been implemented, which automatically adds the generated figures in tasks to the paper dependency list and returns the target path for the tasks targets.
```
target = paper_plot_target('figure_name.pdf')
```

### GitHub action
The pydoit workflow is included in the action, so we have a direct feedback if everything is working as planned.

### tex-macros
In **tex/macros/tex_macros.tex** are macros defined which are static.
This whould be for example the definition of variables in equations, for example
```
\newcommand{\DOH}{\alpha} % degree of hydration
```

### py-macros
There is the possibility to generate tex-macros by a python script.
This allows variable input to be included in the paper.
This could be the computation time, a calculated KPI or the name of a figure included in the paper.
Currently, these macros are defined in **tex/macros/py_macros.yaml**.
It is a simpel dictionary where the key is the tex-command and the value, the tex-code.
The script **tex/macros/py_macros.py** then generates the **tex/macros/py_macros.tex** from this.
The dictionary is also read in the dodo file, so file names can be passed as inputs to functions.

## Comments from Atul:
### Analyze_kpis
This folder has a script to run the snakemake workflow for varrying design variables and saving the KPIs as csv file. Also has script to plot from this .csv file. Careful with the Paths !!
