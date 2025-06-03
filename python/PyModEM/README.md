PyModEM
=======

A small collection of Python utilities to manipulate and plot ModEM Data and
grid files.

# Installation

PyModEM can be installed using `pip`. It's often easier to use a Python virtual
enviorment for installations:

```bash
$ git clone https://github.com/MiCurry/ModEM-Tools.git
$ cd ModEM-Tools/python/PyModEM
$ python -m venv pymodem-venv
$ source pymodem-venv/bin/activate # For Unix/Mac
$ ./pymodem-venv/bin/activate # For Windows
```

Then, install the requirments and PyModEM:

```bash
$ pip install -r requirments.txt
$ pip install -e .
```

[MTPY]: https://github.com/MTgeophysics/mtpy-v2

# Usage

Activate your virtual enviorment, (if you created it above):

```bash
$ source ./pymodem-venv/bin/activate
```

You'll then be able to run the scripts that are found in `ModEM-Tools/python/PyModEM/scirpts`:

```bash
$ modem_data cascad_errfl5.dat list
11.63636
25.6
53.89474
102.4
215.5789
409.6
862.3158
1638.4
4681.143
18724.57
```

You can also create your own files and access PyModEM classes for your own manipulation:

```Python

from PyModEM import ModEMData

data = ModEMData.ModEMData('cascad_errfl5.dat')

print(data.periods)
data.remove_period(25.6)
```

For more examples, view the scripts in `ModEM-Tools/python/PyModEM/scripts`.