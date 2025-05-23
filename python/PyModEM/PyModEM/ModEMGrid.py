import sys
import copy
import math

from typing import List, Tuple
from dataclasses import dataclass

import numpy as np

from PyModEM.ModEMCovariance import ModEMCovariance

def mult_by_one(value, base):
    return value * 1

@dataclass
class ModEMGridData:
    grid_ew : np.ndarray
    grid_ns : np.ndarray
    grid_z : np.ndarray
    rhos : np.ndarray
    origin: Tuple[float, float]
    orientation: 0.0
    mask: np.ndarray

class ModEMGrid:
    def __init__(self, fname=None, data=None):
        if data is not None:
            self._data = data
        else:
            self._data = ModEMGridData(grid_ew=None,
                                       grid_ns=None,
                                       grid_z=None,
                                       rhos=None,
                                       origin=None,
                                       orientation=None,
                                       mask=None)

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

    def create_origin(self):
        self._data.origin = [self.grid_east.min(), self.grid_north.min()]

    def set_grid_north(self, nodes):
        self.grid_north = np.array(
            [
                nodes[0:ii].sum() for ii in range(nodes.size + 1)
            ]
        )

    def parse_header(self, header: str) -> [int, int, str]:
        header_split = header.strip().split()

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


    def read_dim(self, file, n):
        count = 0
        arry = []
        while count != n:
            line = file.readline().strip().split()
            count += len(line)
            arry += line

        arry = np.array(arry)
        return arry.astype(float)

    def read_rhos(self, file, nx, ny, nz):
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

    def load(self, filename: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        with open(filename, 'r') as model_file:
            comment = model_file.readline()
            header = model_file.readline()

            nx, ny, nz, self.resistivity_type = self.parse_header(header)

            if self.resistivity_type == 'LOGE':
                self._convert = math.log
                self._base = math.e
            elif self.resistivity_type == 'LOG10':
                self.convert = math.log
                self._base = math.e
            elif self.resistivity_type == 'LINEAR':
                self.convert = mult_by_one
                self._base = None

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
                                      rhos=self.rhos,
                                      origin=self.origin,
                                      orientation=self.orientation,
                                      mask=None)


    def read_header(self, fname:str):
        with open(fname, 'r') as file:
            self.read_comment(file)
            self.read_dimensions(file)
            self.read_x_dimensions(file)
            self.read_y_dimensions(file)
            self.read_z_dimensions(file)

    def read_comment(self, file):
        self.comment = file.readline()

    def read_dimensions(self, file):
        line = file.readline().strip().split()
        self.nx_read = int(line[0])
        self.ny_read = int(line[1])
        self.nz_read = int(line[2])
        self.ncodes = int(line[3])
        self.resistivity_type = line[4]

    def read_x_dimensions(self, file):
        x_dims = file.readline().strip().split()

        if len(x_dims) != self.nx_read:
            raise ValueError(f"The number of x dimensions was not the same size as specified in line 2 "\
                             f"(got {len(x_dims)} but expected {self.nx})")

        self.x_dims = np.array(x_dims).astype(float)

    def read_y_dimensions(self, file):
        y_dims = file.readline().strip().split()

        if len(y_dims) != self.ny_read:
            raise ValueError(f"The number of y dimensions was not the same size as specified in line 2 "\
                             f"(got {len(y_dims)} but expected {self.nx})")
        
        self.y_dims = np.array(y_dims).astype(float)

    def read_z_dimensions(self, file): 
        z_dims = file.readline().strip().split()

        if len(z_dims) != self.nz_read:
            raise ValueError(f"The number of z dimensions was not the same size as specified in line 2 "\
                             "(got {len(z_dims)} but expected {self.nx})")
        
        self.z_dims = np.array(z_dims).astype(float)

    def read_resistivity_codes(self, file):
        pass

    def read_cell_resistivities(self, file):
        pass

    def read_coordinate_origin(self, file):
        pass

    def read_coordinate_rotation(self, file):
        pass


    def write_grid(self, fname):
        with open(fname, 'w') as file:
            self.write_header_line(file)
            self.write_grid_dimensions_line(file, len(self.grid_ew), len(self.grid_ns), len(self.grid_z))
            self.write_grid_distancing(file, self.grid_ns)
            self.write_grid_distancing(file, self.grid_ew)
            self.write_grid_distancing(file, self.grid_z)
            self.write_resistivity_codes(file)
            self.write_rhos(file, self.rhos)
            self.write_origin(file, self.origin)
            self.write_orientation(file, self.orientation)


    def write_header_line(self, file, comment="# Data writen by PyModEM"):
        file.write(comment + '\n') 

    def write_grid_dimensions_line(self, file, nx: int, ny: int, nz: int, 
                                   nCodes=0,
                                   resitivity_type="LOGE"):

        file.write(f"{ny} {nx} {nz}     {nCodes} {resitivity_type}\n")

    def write_grid_distancing(self, file, distances: np.ndarray):
        nWriten = 0
        for distance in distances:
            if nWriten == 0:
                file.write(f"{abs(distance):>12.3f}")
            else:
                file.write(f" {abs(distance):>12.3f}")

            nWriten += 1

        file.write('\n')

    def write_resistivity_codes(self, file):
        # ModEM does not use resitivity codes
        print("WARNING: Not writing resistivity codes, ModEM does not use resistivity_codes", file=sys.stderr)
        pass

    def write_rhos(self, file, rhos):
        rhos = np.flip(rhos, 0)
        for zz in range(self.grid_z.size):
            file.write("\n")
            for yy in range(self.grid_ew.size):
                for xx in range(self.grid_ns.size):
                    file.write(f"{rhos[xx, yy, zz]:>13.5E}")
                file.write("\n")

    def write_origin(self, file, origin: Tuple[float, float]):
        file.write(f"\n{origin[0]:>16.3f} {origin[1]:16.3f} 0.0\n")

    def write_orientation(self, file, orientation: float):
        file.write(f"{orientation}\n")

    def mask_water(self, water_cond=0.3, esp=0.0005):
        self.mask = ModEMCovariance.mask_water(self.rhos, self.convert(water_cond), esp=esp)

    def write_covariance(self, fname,
                         x_smoothing=0.0,
                         y_smoothing=0.0,
                         z_smoothing=0.0,
                         n_smoothing=0):

        ModEMCovariance.write_covariance(fname,
                                         self,
                                         x_smoothing,
                                         y_smoothing,
                                         z_smoothing,
                                         n_smoothing,
                                         mask=self.mask)

    def convert(self, value):
        return self._convert(value, self._base)


if __name__ == "__main__":
    import argparse
    description = '''ModEM Grid Reader'''

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('rho_filename', help='Filepath to the rho model', type=str)
    args = parser.parse_args()
    
    grid = ModEMGrid(args.rho_filename)
