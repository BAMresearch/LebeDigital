import pathlib
from doit import get_var
from doit.action import CmdAction
from doit.tools import config_changed

def task_paper():
    """compile pdf from latex source"""
    paper = "tex/optimization_paper.tex"
    return {
        "file_dep": [paper],
        "actions": [f"tectonic {paper}"],
        "targets": ["optimization_paper.pdf"],
        "clean": True,
    }
