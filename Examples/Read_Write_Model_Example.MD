Read/Write ModEM Model Example
==============================

We can read and write ModEM Model datafiles by using [readCond_3D][readCond_3D]
and [writeCond_3D][writeCond_3D].

These files can read and write either Weerachai Siripunvaraporn's "0" (WS3D)
format or Mackie's format. By default, ModEM uses the WS3D format, so this
example will explain reading WS3D files.

Before we begin, it's worth nothing that the WS3D format stores data as
resistivity and *not* as conductivity. Thus, read and write conductivity will
convert resistivity to conductivity.

For this example, we will use the Cascadia example that is found in the
[ModEM-Examples][ModEM-Examples] repository. 

[readCond_3D]: ../matlab/ioAscii/readCond_3D.m
[writeCond_3D]: ../matlab/ioAscii/writeCond_3D.m
[ModEM-Examples]: https://github.com/MiCurry/ModEM-Examples

## Reading Conductivity Data

Let us read the `cascad_half_inverse.ws` file from the Cascadia Magnetotelluric
example in the ModEM-Examples repository. We will need to specify the format as 2 to
read in WS3D files, but if we have a file in Mackie format we can pass 1 to format:

```matlab
>> cond = readCond_3D('ModEM-Examples/Magnetotelluric/3D_MT/Cascadia/cascad_half_inverse.ws', 2)

cond = 

  struct with fields:

    paramType: 'LOGE'
            v: [48x46x34 double]
      AirCond: -23.0259
         grid: [1x1 struct]
```

`readCond_3D` returns a struct with four different fields. The `v` field
contains the resistivity values that have been converted into conductivity.

`paramYype` represents the resistivity type that was listed in the model file.
Either: LOGE, LOG10 or LINEAR.

`cond.AirCond` will either read the conductivity found data file, or, for WS3D
files use the default value of paramType(1e-10) (i.e. log10(1e-10), ln(1e-10)
or 1e-10).

`cond.grid` is a structured that defines all information to describe the grid
and it's dimensions. Let's take a look further at `cond.grid`.

### Conductivity Grid Structure (`cond.grid`)

```matlab

>> cond.grid
ans = 

  struct with fields:

          dx: [48x1 double]
          dy: [46x1 double]
          dz: [34x1 double]
          Nx: 48
          Ny: 46
     NzEarth: 34
       NzAir: 6
      origin: [3x1 double]
    rotation: 0
       units: 'km'
```

`dx`,  `dy` and `dz` contain the `distances` of the grid cells in their
respective directions in the units defined in `grid.units` (in kilometers/km).

As you can see, `Nx`, `Ny` and `NzEarth` contain the number of cells in the x,
y, and z direction. Likewise, `NzAir` contains the number of air layers defined
in the file. WS3D files do not have the ability to define air layers, thus it
defaults to 6.

The origin is important as that it represents the origin of the relationship to
the data file. For more information, see the 
[ModEM Users Guide][modem-users-guide] section 4.2.

[modem-users-guide]: https://github.com/MiCurry/ModEM-Model/blob/main/doc/userguide/ModEM_UserGuide.pdf

## Writing Conductivity Data

Continuing from the example before, we can easily write out conductivity values
by using `writeCond_3D`. To write a WS3D file we can specify it by using passing
2 to format in `writeCond_3D`.

```matlab
>>> status = writeCond_3D('example.cascade.ws3d.rho', cond, 2)

status =

     0
```

Status will contain the `status` output argument returned from MatLab's `fclose`
(i.e. 0 if successful, -1 on error).

