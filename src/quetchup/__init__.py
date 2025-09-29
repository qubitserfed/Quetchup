from quetchup.Extension import *
import IPython.display

def to_latex(state : State, ipynb=True) -> str:
    if ipynb:
        return IPython.display.Math(state.__latex_rep())
    else:
        return state.__latex_rep()

def tensor(a, b):
    if isinstance(a, State) and isinstance(b, State):
        return state_tensor(a, b)
    elif isinstance(a, Map) and isinstance(b, Map):
        return map_tensor(a, b)
    else:
        raise ValueError("Invalid input types")
