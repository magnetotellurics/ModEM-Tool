Read/Write ModEM Data MatLab Example
====================================

This example will walk you through reading and writing ModEM 'list'/ASCII
datafiles using the ModEM MatLab tools.

For this example, we will read and write a few different data types. This
example will use examples found in the [ModEM-Examples][ModEM-Examples].

[ModEM-Examples]: https://github.com/MiCurry/ModEM-Examples

All of the functions listed in this example contain usage documentation, which
you can see at the top of the file, or by typing `help function` in MatLab.

## Reading Impedance Data types using readZ_3D

To read impedance data types, we can use [readZ_3D][readZ_3D]:

```matlab
>> [data, header, units, isign, origin, info] = readZ_3D('./ModEM-Examples/Magnetotelluric/3d_MT/Cascadia/cascad_errfl3.dat');
```

readZ_3D can read in the following datatypes:

* Full_Impedance
* Off_Diagonal_Impedance
* Full_Vertical_Components
* Off_Diagonal_Rho_Phase - Will be converted into Full_Impedance - Use
 [readApres_3D][readApres_3D] to read in Off_Diagonal_Rho_Phase without converting

As noted above, you can also read in Off_Diagonal_Rho_Phase data, but it will be
converted into Full_Impedance. If you want to read in Off_Diagonal_Rho_Phase
without converting you can use [readApres_3D][readApres_3D].

readZ_3D returns a number of output arguments, the most which are self
explanatory; however, data and info are data structures which are worth
explaining.

### allData data structure

```matlab
>> data

data =

  1x10 cell array

    {1x1 struct}    {1x1 struct}    {1x1 struct}    {1x1 struct}    {1x1 struct}    {1x1 struct}    {1x1 struct}    {1x1 struct}    {1x1 struct}    {1x1 struct}
```

The `allData` output argument returned from [readZ_3D][readZ_3D] is MatLab cell
array that contains one cell per period for all data types that were found in
the data file.

We can take a look at the data structure of each cell by accessing one cell:

```matlab
>> data{1}

ans = 

  struct with fields:

                 T: 11.6364
             Cmplx: 1
    signConvention: 1
             nComp: 12
          nanvalue: NaN
           siteLoc: [109x3 double]
          siteChar: {1x109 cell}
                 Z: [109x6 double]
              Zerr: [109x6 double]
            origin: [45.2760 -119.6340 0]
            orient: 0
               lat: [109x1 double]
               lon: [109x1 double]
              type: ["Full_Impedance"    "Full_Vertical_Components"]
          compChar: [6x3 char]
             units: '[]'
```

Here we can see that two data types have been read, `Full_Impedance` followed
by `Full_Vertical_Components`. The data values for these two data types exist
in side of `data{1}.Z` in the order that they are listed in the `type` field.

Likewise, the `compChar` contains the order of the data components in the Z and
Zerr field:

```matlab
>> data{1}.compChar

ans =

  6x3 char array

    'ZXX'
    'ZXY'
    'ZYX'
    'ZYY'
    'TX '
```

> *Note:* readZ_3D will always return `Full_Impedance` followed by
> `Full_Vertical_Components` in the `allData` cell arrays even if they appear in
> the reverse order in the data file.
>
> This is to ensure that older scripts can use these two datatypes, which often expect
> `Full_Impedance` followed by `Full_Vertical_Components`.

### Info Data Structure

The `info` returned by readZ_3D is a cell array which each array element being
associated with each data type found in the file. In the example above, we can
see that there `info` is a 1x2 cell array:

```matlab
>> info

info =

  1x2 cell array

    {1x1 struct}    {1x1 struct}
```

If we access the first cell, we can see it contains information for Full_Impedance:

```matlab
>> info{1}

ans = 

  struct with fields:

     data: [109x10x4 double]
      err: [109x10x4 double]
     type: 'Full_Impedance'
    units: '[mV/km]/[nT]'
      lat: [41.8760 48.9350 47.0440 46.3530 45.8890 45.7670 45.3990 44.6460 44.7140 44.5580 44.0520 43.9170 43.4630 43.4360 42.7750 42.7320 42.0790 42.1530 44.6680 45.9330 45.7940 ... ] (1x109 double)
      lon: [-124.1900 -116.4440 -116.3460 -116.2130 -116.1570 -115.2230 -116.2790 -116.6920 -116.0100 -114.8470 -116.8300 -115.9620 -116.7900 -115.9370 -116.8690 -115.8680 -116.4770 ... ] (1x109 double)
      loc: [109x3 double]
     code: [109x5 char]
      per: [10x1 double]
    ncomp: 8
     comp: [4x3 char]
```

`info{2}` contains information for Full_Vertical_Components. But as you can see
info structure contains all information for the associated data type.

[readZ_3D]: ../matlab/ioAscii/readZ_3D.m

## Reading Phase data types using readApres_3D

To read Phase_Tensor and Off_Diagonal_Rho_Phase data types, we need to use 
[readApres_3D][readApres_3D]. 

```matlab
function [allData, header, units, isign, origin, info] = readApres_3D(cfile, newunits, onetype)
```

`readApres_3D` can read non-impedance data types:

* Off_Diagonal_Rho_Phase
* Phase_Tensor

As you can see, the input arguments and output arguments are roughly the same as
readZ_3D. `allData` returned by `readApres_3D` is the same as what is described
above. i.e. a cell array with one cell per period.

Likewise, `info` is a cell array that with one cell per data type read inside the file.

[readApres_3D]: ../matlab/ioAscii/readApres_3D.m

## Reading data files that contain both impedance and apparent resistivity data types

If you have both impedance and apparent resistivity data types inside a single
file, you can read them both in by using both `readZ_3D` and `readApres_3D`
and using the `onetype` argument.

```matlab
[imp_data, imp_header, imp_units, imp_sign, imp_origin, imp_info] = readZ_3D('full_imp_n_phase.dat', '', 'Full_Impedance');
[phase_data, phase_header, phase_units, phase_sign, phase_origin, phase_info] = readApres_3D('full_imp_n_phase.dat', '', 'Phase_Tensor');
```

# Write Data with writeZ_3D  

Any datatype can be written out with [writeZ_3D][writeZ_3D]:

```matlab
function [status] = writeZ_3D(cfile, allData, header, units, isign, convert_to_apres)
```

You can pass the `allData` that is returned from `readZ_3D` to `writeZ_3D` and
by specifying the filename in `cfile`:

```matlab
>> [data, header, units, isign, origin, info] = readZ_3D('~/Projects/ModEM-Examples/Magnetotelluric/3d_MT/Cascadia/cascad_errfl3.dat');
>> [status] = writeZ_3D('example.dat', data);
```

writeZ_3D will write out both Full_Impedance and Full_Vertical found in
`cascad_errfl3.dat`.

[writeZ_3D]: ../matlab/ioAscii/writeZ_3D.m
