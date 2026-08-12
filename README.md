# ModEM-Tools

The ModEM-Tools repository contains a collection of MATLAB and Python tools for
manipulating ModEM data and model files.

**All of the MatLab scripts and classes in this repository have been created and developed by Dr. Anna Kelbert.**

This repository has the following folders:

**matlab-tools**

Still under development, these `matlab-tools` directory contains classes for
handling and manipulating electromagnetic data, reading magnetotelluric (MT)
transfer functions (reading SPUD bundles) and classes for geospatial
transformations of model parameterizations. 

**modem-io**

A collection of MATLAB and Python code for reading, writing and plotting ModEM
data types. The MATLAB code has been more widely used and is more recommended
then the Python code.

**scripts** 

Contains scripts that are not longer used, but may act as examples.

**utilities**

The utilities folder contains miscellaneous MATLAB classes such as a GUI for
viewing 2D and 3D model files, tools for processing EDI files and others.

## Note on these scripts/tools

Some of these tools **may** be incomplete, contain errors or bugs, or may
not work at all. We are providing these tools as a starting point or reference
for derivative or future work.

However, we have done work to ensure that some basic functionality is working.
This functionality can be found in the [MATLAB Examples](#matlab-examples)
section of this README.

## Related Repositories

* [ModEM-Model][ModEM-Model] - The ModEM Model itself
* [ModEM-Examples][ModEM-Examples] - A collection of MT and CSEM examples for ModEM

[ModEM-Model]: https://github.com/magnetotellurics/ModEM
[ModEM-Examples]: https://github.com/magnetotellurics/ModEM-Examples


## MATLAB Code

> **NOTE:** Some MATLAB code requires the MATLAB [Mapping Toolbox][Mapping ToolBox].

The collection of MATLAB code in this repository contain a number of functions
and classes that have been used for many years by Anna Kelbert and have a
plethora of functionality. It is the best resource for tools for manipulating
ModEM data, model files or performing other actions.

[Mapping ToolBox]: https://www.mathworks.com/products/mapping.html

### MATLAB Examples

We have provided a few examples using MATLAB code which should help anyone
getting started with ModEM datatypes or model files. 

* [Read/Write ModEM Data Using readZ_3D and writeZ_3D][rw_data_example] - A
basic example using the low-level `readZ_3D` and `writeZ_3D` functions to read
and write  a data file
* [Read/Write ModEM Model Example][rw_model_example] - A basic example using the
low-level to read and write a model file
* [Using MATLAB classes to create a synthetic example][matlab_classes_example] - An example
that uses the higher level MATLAB classes to create synthetic data and a
synthetic model as well as an example using the classes to plot a model.

[rw_data_example]: /modem-io/Examples/Read_Write_Data_Example.MD
[rw_model_example]: /modem-io/Examples/Read_Write_Model_Example.MD
[matlab_classes_example]: /modem-io/Examples/Creating_Synthetic_Model_and_Data.md

### modem-io MATLAB classes

The MATLAB classes provide the greatest functionality and are very flexible.
They will most likely be the best option to meet your needs. 

Some of these classes are used 

* [xygrid][xygrid] and [llgrid][llgrid] - Classes for manipulating (or creating) a grid
* [xymodel][xymodel] and [llmodel][llmodel] - Class for manipulating a gird and it's conductivity field
* [mtdata][mtdata] - Class for manipulating a ModEM data file (contains a collection of `mtperiods`)
* [mtperiod][mtperiod] - Class for manipulating and defining a period, including station data.
* [mttf][mttf] - Read in either Z, XML, EDI or BIRRP files

If you want to use these classes to read in an existing file, make sure you check out their
static methods. For instance:

```matlab
>>> data = mtdata.read('datafile.dat');
>>> model = xymodel.read('example.rho');
```

## Python Code (PyModEM)

The Python Code is a newer tool with significantly less functionality then the
MATLAB code. The MATLAB code should be used in favor of it, but we are providing
it here as reference and as a starting point for future tools.

Currently it has some basic functionality that allow for the creation of
synthetic data, model files and covariance files as well as some terse 
plotting functionality.

Some of the scripts require the [MtPy][MtPy] Python tool. MtPy is a great 
resource and we recommend using it as well.

For more information on the PyModEM please see it's [README][pymod-readme].

[MtPy]: https://github.com/MTgeophysics/mtpy-v2
[pymod-readme]: /modem-io/PyModEM/README.md
