import glob
from runipy.notebook_runner import NotebookRunner
from IPython.nbformat.current import read

from easydev.console import purple
from easydev import Progress

def check_ipython_notebook():


    notebooks = glob.glob("*ipynb")
    N = len(notebooks)

    pb = Progress(N)
    for i,filename in enumerate(notebooks):
        print(purple(filename))
        notebook = read(open(filename), 'json')
        r = NotebookRunner(notebook)
        r.run_notebook()
        pb.animate(i+1)



if __name__ == "__main__":
    check_ipython_notebook()
