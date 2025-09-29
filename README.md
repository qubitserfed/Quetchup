![Quetchup Logo](logo.png)

# Quetchup
A python package for **phase accurate** simulation of the Clifford category. Can decide equality (and more!) for all linear maps obtained by composing or tensoring the following:
* Clifford Unitaries
* Projections on $| 0 \rangle$
* $|0 \rangle$ state preparation

## Installation
This is a python package, so just clone the repo and then `pip install [the repo location on your computer]`.

## Usage
```python
from quetchup import *
x = x_map(), z = z_map(), swap = swap_map()
print(swap * tensor(x, z) * swap == tensor(z, x))
```

For more examples, see the [tutorial](examples/tutorial.ipynb) notebook.
