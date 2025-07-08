import sys
import copy
import math

from typing import Tuple, TextIO
from dataclasses import dataclass

import numpy as np

from PyModEM.ModEMCovariance import ModEMCovariance

def mult_by_one(value, base):
    return value * 1

@dataclass
class ModEMGridData:
    ''' ModEMGridData - Holds information on ModEM Grid/rho file

    This class should be used in conjuction with ModEMGrid class.

    Parameters
    ----------
    grid_ew : np.ndarray
        Grid spacing in the east, west direction
    grid_ns : np.ndarray
        Grid spacing in the north, south direction
    grid_z : np.ndarray
        Grid spacing in the z direction
    rhos : np.ndarray
        Cell resistivities
    origin : Tuple[float, float]
        Coordinate origin
    orientation : float
        coordinate rotation (not used in ModEM)
    mask : np.ndarray
        Mask of cell resistivies to use to create covariance mask
    '''
    grid_ew : np.ndarray
    grid_ns : np.ndarray
    grid_z : np.ndarray
    resistivity_type : str
    rhos : np.ndarray
    origin: Tuple[float, float]
    orientation: float
    mask: np.ndarray

class ModEMGrid:
    def __init__(self, fname=None, data:ModEMGridData=None):
        ''' ModEMGrid - A class used to create, read, edit, alter a ModEM Rho/Sigma/Model file

        If fname is present, instantiate the class with the model file found at fname.

        If data is present, then instantiate the class with the data inside data.

        If fname, or data are both not present, then instantiate a class with all ModEMGridData
        fields to be None.

        Parameters
        ----------
        fname : str
            Filename to read a ModEM grid from
        data : ModEMGridData
            ModEMGridData object to use to instantiate this class

        '''
        if data is not None:
            self._data = data
            self._determine_conversion_function(data.resistivity_type)

        else:
            self._data = ModEMGridData(grid_ew=None,
                                       grid_ns=None,
                                       grid_z=None,
                                       resistivity_type='LOGE',
                                       rhos=None,
                                       origin=None,
                                       orientation=None,
                                       mask=None)
            self._determine_conversion_function(self._data.resistivity_type)

        if fname is not None:
            self.fname = fname
            self.load(fname)
            data=self._data

    def __copy__(self):
        return ModEMGrid(data=copy.copy(self._data))

    def __deepcopy__(self, memo):
        return ModEMGrid(data=copy.deepcopy(self._data, memo))

    @property
    def origin(self) -> Tuple[float]:
        return self._data.origin

    @origin.setter
    def origin(self, value):
        self._data.origin = value

    @property
    def orientation(self) -> float:
        return self._data.orientation

    @orientation.setter
    def orientation(self, value) -> float:
        self._data.orientation= value

    @property
    def nx(self) -> int:
        return len(self._data.grid_ew)

    @property
    def neast(self) -> int:
        return self.nx()

    @property
    def ny(self) -> int:
        return len(self._data.grid_ns)

    @property
    def nnorth(self) -> int:
        return self.ny()

    @property
    def nz(self) -> int:
        return len(self._data.grid_z)

    @property
    def ndepth(self) -> int:
        return self.nz()

    @property
    def grid_ew(self) -> np.ndarray:
        return self._data.grid_ew

    @property
    def grid_ns(self) -> np.ndarray:
        return self._data.grid_ns

    @property
    def grid_z(self) -> np.ndarray:
        return self._data.grid_z

    @property
    def rhos(self) -> np.ndarray:
        return self._data.rhos

    @rhos.setter
    def rhos(self, value):
        self._data.rhos = value

    @property
    def mask(self) -> np.ndarray:
        return self._data.mask

    @mask.setter
    def mask(self, value):
        self._data.mask = value

    def _determine_conversion_function(self, res_type):
        if res_type == 'LOGE':
            self._log_function = math.log
            self._base = math.e
        elif res_type == 'LOG10':
            self._log_function = math.log
            self._base = math.e
        elif res_type == 'LINEAR':
            self._log_function = mult_by_one
            self._base = None
        else:
            raise ValueError(f"Did not understand resistivity type {res_type} choose: [ 'LOGE', 'LOG10', 'LINEAR' ]")

    @property
    def resistivity_type(self) -> str:
        return self._data.resistivity_type

    @resistivity_type.setter
    def resistivity_type(self, value):
        self._determine_conversion_function(value)
        self._data.resistivity_type = value

    def create_origin(self):
        self._data.origin = [self.grid_east.min(), self.grid_north.min()]

    def set_grid_north(self, nodes):
        self.grid_north = np.array(
            [
                nodes[0:ii].sum() for ii in range(nodes.size + 1)
            ]
        )

    def parse_dimension_line(self, dimension_line: str) -> Tuple[int, int, int, str]:
        ''' parse_dimension_line(self, header) - Parse the Dimension line
        string. This method will read grid dimensions (nx, ny, nz),  number of
        resistivity codes (not used) and resistivity type (LINEAR, LOG10, LOGE))

        Should be in the format of:

        ```
        nx ny nz  ncodes res_type
        ```

        Parameters
        ----------
        dimension_line : str
            String that contains the nx, ny, nz, ncodes and resistivity type
        '''
        header_split = dimension_line.strip().split()

        if len(header_split) != 5:
            raise ValueError(f"{self.fname} does not appear to be a 3D Resistivity model in Randi Mackie's format.")

        try:
            nx = int(header_split[0])
        except ValueError as error:
            raise ValueError(f"nx '{header_split[0]}' needs to be an integer value") from error

        try:
            ny = int(header_split[1])
        except ValueError as error:
            raise ValueError(f"ny '{header_split[1]}' needs to be an integer value") from error

        try:
            nz = int(header_split[2])
        except ValueError as error:
            raise ValueError(f"nz '{header_split[2]}' needs to be an integer value") from error

        try:
            resistivity_codes = int(header_split[3])
        except ValueError as error:
            print(f"INFO: resistivty code was not an integer (but is ignored by modem) skipping..")

        try:
            resistivity_type = str(header_split[4])
        except ValueError as error:
            raise ValueError(f"Resitivity type '{header_split[4]}' needs to be an string value") from error

        if resistivity_type != 'LOGE' and resistivity_type != 'LINEAR' and resistivity_type != 'LOG10':
            raise ValueError(f"Error: Unepexcted resistivity type: {resistivity_type}")

        return nx, ny, nz, resistivity_type 


    def read_dim(self, file: TextIO, n: int) -> np.ndarray:
        ''' read_dim(self, file, n) - Read model file grid spacing/dimensions

        Parameters
        ----------
        file : TextIO
            Open file to read from
        n : int
            Size/number of dimensions/cells in this direction

        Returns
        -------
        out : np.ndarray
            Grid dimension spacing in this direction
        '''
        count = 0
        arry = []
        while count != n:
            line = file.readline().strip().split()
            count += len(line)
            arry += line

        arry = np.array(arry)
        return arry.astype(float)

    def read_rhos(self, file: TextIO, nx: int, ny: int, nz: int) -> np.ndarray:
        ''' read_rhos(self, file, nx, ny, nz) - Read the rhos/resistivites from a file

        Parameters
        ----------
        file : TextIO
            Open file to read from
        nx : int
            Number of cells in the x direction
        ny : int
            Number of cells in the y direction
        nz : int
            Number of cells in the z direction

        Return
        ------
        out : np.ndarray
            Numpy array holding the resistivities of the file in shape [nx, ny, nz]

        '''
        count = 0
        count_e = 0
        rhos = np.zeros([nx, ny, nz])

        while count < nz:
            spline = file.readline().strip().split()

            if len(spline) == 0:
                count += 1
                count_e = 0
            elif (len(spline) == 3) & (count == nz - 1):
                count += 1 
                count_e = 0
            else:
                north_line = np.array([float(nres) for nres in spline])

                rhos[:, count_e, count] = north_line[::-1]
                count_e += 1

            if len(spline) == 3:
                self.origin = spline

        return rhos

    def read_origin(self, file):
        return np.array(self.origin).astype(float)

    def read_orientation(self, file):
        return float(file.readline().strip().split()[0])

    def load(self, filename: str):
        ''' ModEMGrid.load(filename) - Read a model file into this class

        Parameters
        ----------
        filename : str
            File name to open and read
        '''
        with open(filename, 'r') as model_file:
            comment = model_file.readline()
            header = model_file.readline()

            nx, ny, nz, resistivity_type = self.parse_dimension_line(header)

            self._determine_conversion_function(resistivity_type)

            self.x = self.read_dim(model_file, nx)
            self.y = self.read_dim(model_file, ny)
            self.z = self.read_dim(model_file, nz)
            
            tmp = model_file.readline()

            self.rhos = self.read_rhos(model_file, nx, ny, nz)
            self.origin = self.read_origin(model_file)
            self.orientation = self.read_orientation(model_file)

            self._data = ModEMGridData(grid_ew=self.y,
                                      grid_ns=self.x,
                                      grid_z=self.z,
                                      resistivity_type=resistivity_type,
                                      rhos=self.rhos,
                                      origin=self.origin,
                                      orientation=self.orientation,
                                      mask=None)

    def write_grid(self, fname : str):
        ''' write_grid(self, fname) - Write out the current object in the WS/Modem format 
        
        Parameters
        ----------
        fname : str
            Write out the current object as a WS/Modem Format to fname

        '''
        with open(fname, 'w') as file:
            self._write_header_line(file)
            self._write_grid_dimensions_line(file, len(self.grid_ew), len(self.grid_ns), len(self.grid_z), resitivity_type=self.resistivity_type)
            self._write_grid_distancing(file, self.grid_ns)
            self._write_grid_distancing(file, self.grid_ew)
            self._write_grid_distancing(file, self.grid_z)
            self._write_resistivity_codes(file)
            self._write_rhos(file, self.rhos)
            self._write_origin(file, self.origin)
            self._write_orientation(file, self.orientation)


    def _write_header_line(self, file, comment="# Data writen by PyModEM"):
        file.write(comment + '\n') 

    def _write_grid_dimensions_line(self, file, nx: int, ny: int, nz: int,
                                   nCodes=0,
                                   resitivity_type="LOGE"):

        file.write(f"{ny} {nx} {nz}     {nCodes} {resitivity_type}\n")

    def _write_grid_distancing(self, file, distances: np.ndarray):
        nWriten = 0
        for distance in distances:
            if nWriten == 0:
                file.write(f"{abs(distance):>12.3f}")
            else:
                file.write(f" {abs(distance):>12.3f}")

            nWriten += 1

        file.write('\n')

    def _write_resistivity_codes(self, file):
        # ModEM does not use resitivity codes
        print("WARNING: Not writing resistivity codes, ModEM does not use resistivity_codes", file=sys.stderr)
        pass

    def _write_rhos(self, file, rhos):
        rhos = np.flip(rhos, 0)
        for zz in range(self.grid_z.size):
            file.write("\n")
            for yy in range(self.grid_ew.size):
                for xx in range(self.grid_ns.size):
                    file.write(f"{rhos[xx, yy, zz]:>13.5E}")
                file.write("\n")

    def _write_origin(self, file, origin: Tuple[float, float]):
        file.write(f"\n{origin[0]:>16.3f} {origin[1]:16.3f} 0.0\n")

    def _write_orientation(self, file, orientation: float):
        file.write(f"{orientation}\n")

    def mask_water(self, water_cond=0.3, esp=0.0005) -> np.ndarray:
        ''' mask_water(self, water_cond=0.3, esp=0.0005) - Mask water where rho
        value equals water_cond +/- esp

        Create a ModEM covariance mask where rho values that are:

            water_cond - esp <= rho <= water_cond + esp

        All values that meet the above condition will be set to 9 in the mask
        and all other values will be set to 1.

        Values of water should be in ohm-meters and will be converted using the
        resistivity type that is set in self.resistivity_type to match the
        values of rho.

        Parameters
        ----------
        water_cond : float

        water_cond : float
            Value to use f

        esp : float
            epsiolon value to use to help with matching values 'near' water_cond:

            water_cond - esp >= rho >= water_cond + esp

        Returns
        -------
        Array of masked values
        '''
        self.mask = ModEMCovariance.mask_water(self.rhos, self.convert(water_cond), esp=esp)
        return self.mask

    def write_covariance(self, fname: str,
                         x_smoothing=0.0,
                         y_smoothing=0.0,
                         z_smoothing=0.0,
                         n_smoothing=0):
        ''' write_covariance(fname, x_smoothing, y_smoothing, z_smoothing, n_smoothing) - Create covariance file

        Will create a ModEM covariance file with provided x/y/z smoothing. If
        self.mask, if present, the mask will be written out.

        Currently does not provide any ability to write out covariance rules.

        Parameters
        ----------
        fname : str
            file to write
        x_smoothing : float
            Smoothing to apply in the x direction (sets all values of x smoothing to x_smoothing)
        y_smoothing : float
            Smoothing to apply in the y direction (sets all values of y smoothing to y_smoothing)
        z_smoothing : float
            Smoothing to apply in the z direction (sets all values of z smoothing to z_smoothing)
        n_smoothing : int
            The number of times to apply smoothing
        '''

        ModEMCovariance._write_covariance(fname,
                                         self,
                                         x_smoothing,
                                         y_smoothing,
                                         z_smoothing,
                                         n_smoothing,
                                         mask=self.mask)

    def convert(self, value: float) -> float:
        ''' convert(self, value) - Convert value, to the resistivity type
        specified for this Grid

        Wrapper function to ease conversion of values into the resitivity type
        that is specififed for this grid.

        Parameters
        ----------
        value : float
            Value to convert
        '''
        return self._log_function(value, self._base)
