"""Run the MR practical's calculations without Jupyter."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def display(value):
    print(value.to_string() if hasattr(value, 'to_string') else value)


if __name__ == '__main__':
    notebook = Path(__file__).with_name('MR_practical.ipynb')
    cells = json.loads(notebook.read_text())['cells']
    scope = {'display': display}
    # Keep plots in the HTML handout; this entry point prints the calculations.
    plt.show = lambda: None
    for cell in cells:
        source = ''.join(cell['source'])
        if cell['cell_type'] == 'markdown' and source.startswith('## Q'):
            print('\n' + source.splitlines()[0].lstrip('# '))
        elif cell['cell_type'] == 'code':
            exec(compile(source, str(notebook), 'exec'), scope)
    plt.close('all')
