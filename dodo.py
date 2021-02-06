import glob

def task_website():
    pages = glob.glob("doc/*.md")
    return {
        "file_dep": pages + ["conf.py", "index.rst"],
        "actions": ["sphinx-build . .build"],
        "verbosity": 2,
    }
