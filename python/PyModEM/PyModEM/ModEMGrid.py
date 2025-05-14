from dataclasses import dataclass
from typing import List, Tuple
import numpy as np

@dataclass
class ModEMGridData:
    grid_ew : np.ndarray
    grid_ns : np.ndarray
    grid_z : np.ndarray
    rhos : np.ndarray
    origin: Tuple[float, float]
    orientation: 0.0


class ModEMGrid:
    def __init__(self, fname=None, data=None):
        if fname is not None:
            self.fname = fname
            self.read_header(fname)
        
        if data is not None:
            self._data = data
        else:
            self._data = None

    def create_origin(self):
        self._data.origin = [self.grid_east.min(), self.grid_north.min()]
        print("Origin:", self._data.origin)


    def set_grid_north(self, nodes):
        self.grid_north = np.array(
            [
                nodes[0:ii].sum() for ii in range(nodes.size + 1)
            ]
        )

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
    def origin(self) -> Tuple[float]:
        return self._data.origin

    @property
    def orientation(self) -> float:
        return self._data.orientation

    def write_grid(self, fname):
        with open(fname, 'w') as file:
            self.write_header_line(file)
            self.write_grid_dimensions_line(file, len(self.grid_ew), len(self.grid_ns), len(self.grid_z))
            self.write_grid_distancing(file, self.grid_ew)
            self.write_grid_distancing(file, self.grid_ns)
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
        pass

    def write_rhos(self, file, rhos):
        print(self.grid_z.shape[0], self.grid_z.size, self.rhos.shape)
        for zz in range(self.grid_z.size):
            file.write("\n")
            for ee in range(self.grid_ew.size):
                for nn in range(self.grid_ns.size):
                    file.write(f"{rhos[nn, ee, zz]:>13.5E}")
                file.write("\n")

    def write_origin(self, file, origin: Tuple[float, float]):
        file.write(f"\n{origin[0]:>16.3f} {origin[1]:16.3f} 0.0\n")

    def write_orientation(self, file, orientation: float):
        file.write(f"{orientation}\n")

    def write_cov():
        pass


if __name__ == "__main__":
    import argparse
    description = '''ModEM Grid Reader'''

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('rho_filename', help='Filepath to the rho model', type=str)
    args = parser.parse_args()
    
    grid = ModEMGrid(args.rho_filename)
