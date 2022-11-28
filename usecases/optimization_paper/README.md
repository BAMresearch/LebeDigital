## Overview
This directory contains to code to automatically generate the planned paper.
This includes the scripts to generate the figures and at least the option to re-run the computations if needed.

## Structure
The current proposal is the write the paper in the **tex** folder and split the document in sections, which are located in **tex/tex_sections**.
Modules like **workflow_graph** can be added as directory to the main directory and the included in the pydoit workflow.
Figures generated for the paper are located in **figures**.

## Features

### pydoit
There is a **dodo.py** file which includes all steps required to generate the paper.
It has the option for different modes to be implemented.
As example the default is currently
```
doit mode=CI
```
If there are expensive computations the should be skipped, or settings changed for integration testing puroses, this can be done using this option.
### Github action
The pydoit workflow is included in the actions so we have a direct feedback if everything is working as planned.

### tex-macros
In **tex/macros/tex_macros.tex** are macros defined which are static.
This whould be for example the definition of variables in equations, for example
```
\newcommand{\DOH}{\alpha} % degree of hydration
```

### py-macros
There is the possiblity to generate tex-macros by a python script.
This allows variable input to be included in the paper.
This could be the computation time, a calculated KPI or the name of a figure included in the paper.
Currently these macros are defined in **tex/macros/py_macros.yaml**.
It is a simpel dictionary where the key is the tex-command and the value, the tex-code.
The script **tex/macros/py_macros.py** then generates the **tex/macros/py_macros.tex** from this.
The dictionary is also read in the dodo file, so file names can be passed as inputs to functions.
