import sys
import numpy as np
from typing import TextIO

header = "+-----------------------------------------------------------------------------+\n"\
   "| This file defines model covariance for a recursive autoregression scheme.   |\n"\
   "| The model space may be divided into distinct areas using integer masks.     |\n"\
   "| Mask 0 is reserved for air; mask 9 is reserved for ocean. Smoothing between |\n"\
   "| air, ocean and the rest of the model is turned off automatically. You can   |\n"\
   "| also define exceptions to override smoothing between any two model areas.   |\n"\
   "| To turn off smoothing set it to zero. This header is 16 lines long.         |\n"\
   "| 1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)                   |\n"\
   "| 2. Smoothing in the X direction (NzEarth real values)                       |\n"\
   "| 3. Smoothing in the Y direction (NzEarth real values)                       |\n"\
   "| 4. Vertical smoothing (1 real value)                                        |\n"\
   "| 5. Number of times the smoothing should be applied (1 integer >= 0)         |\n"\
   "| 6. Number of exceptions (1 integer >= 0)                                    |\n"\
   "| 7. Exceptions in the form e.g. 2 3 0. (to turn off smoothing between 3 & 4) |\n"\
   "| 8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.|\n"\
   "+-----------------------------------------------------------------------------+\n"

class ModEMCovariance:
    ''' ModEMCovariance - A Static Class to help with creating a ModEM Covariance file.

    At the present time, this class will need to be used with the ModEMData class.
    '''

    AIR_MASK=0
    NO_MASK=1
    WATER_MASK = 9

    def _write_header(file: TextIO ):
        ''' Write the fixed-length header for a ModEM covariance file to an open file

        Parameters
        -----------
        file : TextIO
            Open file or TextIO to write to.
        '''
        file.write(header)
        file.write('\n')

    def _write_grid_dims(file: TextIO,
                        nx: int,
                        ny: int,
                        nz: int):
        ''' Write the grid dimension to an open file (line 1)

        Parameters
        -----------
            file : TextIO
                Open file or TextIO object to write to
            nx, ny, nz : Int
                Grid dimension in nx, ny and nz directions
        '''
        file.write(f"{ny} {nx} {nz}\n")
        file.write("\n")

    def _write_smoothing(file, n, smoothing: float):
        ''' Write the grid dimension to an open file (line 1)

        Parameters
        -----------
            file : TextIO
                Open file or TextIO object to write to
            nx, ny, nz : Int
                Grid dimension in nx, ny and nz directions
        '''
        for i in range(0, n, 1):
            file.write(f"{smoothing} ")
        file.write('\n')

    def _write_modem_smoothing(file,
                              grid,
                              x_smoothing: float,
                              y_smoothing: float,
                              z_smoothing: float,
                              n_smoothing):
        ModEMCovariance._write_smoothing(file, grid.nz, x_smoothing)
        ModEMCovariance._write_smoothing(file, grid.nz, y_smoothing)
        ModEMCovariance._write_smoothing(file, 1, z_smoothing)
        file.write('\n')
        file.write(str(n_smoothing)+'\n')
        file.write('\n')

    def _write_modem_exception_rules(file, n_rules):
        file.write(str(n_rules)+'\n')

    def _write_mask(file, grid, mask):
        for zz in range(grid.grid_z.size):
            file.write(f"{zz+1} {zz+1}\n")
            for xx in range(grid.grid_ns.size):
                for yy in range(grid.grid_ew.size):
                    file.write(f"{mask[xx, yy, zz]} ")
                file.write("\n")

    def write_covariance(fname: str,
                              grid,
                              x_soomthing: float,
                              y_smoothing: float,
                              z_smoothing: float,
                              n_smoothing: int,
                              mask : np.ndarry = None):
        ''' write_covariance(fname, grid, x_smoothing, y_smoothing, z_smoothing, n_smoothing, mask)

        Open and write a ModEM covariance file in fname using grid and the above
        arguments. If mask is present, write out the mask.

        Note: This method and class currently do not have the ability to write ModEM covraiance
        exception rules. You'll need to write them by hand and you might be able to use 'mask_value'
        below, but I have not yet done thourugh tested added multiple masks.

        Parameters
        ----------
        fname : str
            Filename to use to write the covariance file too
        grid : ModEMGrid
            The ModEMGrid you want to write the covariance for.
        x_smoothing : float
            x smoothing to use for all x's in the x direction
        y_smoothing : float
            y smoothing to use for all y's in the y direction
        z_smoothing : float
            z smoothing to use for all z's in the z direction
        n_smoothing : in
            Number of times to apply smoothing
        mask : np.ndarray
            Array of size [grid.nx, grid.ny, grid.nz] to use for writing out the water/land/rule mask
        '''

        with open(fname, 'w') as file:
            ModEMCovariance._write_header(file)

            """ TODO: We should not need to pass in 'grid' as an argument to this function.
            Instead we should pass in a Tuple[nx, ny, nz], but I'm too lazy to change it at the moment.
            """
            ModEMCovariance._write_grid_dims(file, grid.nx, grid.ny, grid.nz)
            ModEMCovariance._write_modem_smoothing(file, grid, x_soomthing, y_smoothing, z_smoothing, n_smoothing)
            ModEMCovariance._write_modem_exception_rules(file, 0)

            if mask is not None:
                ModEMCovariance._write_mask(file, grid, mask)

    def mask_value(rhos, value, mask_value, esp=0.00005, mask=None):
        if mask is None:
            mask = np.full_like(rhos, 1, dtype=int)

        indices = np.where((rhos >= value - esp) & (rhos <= value + esp))

        values_masked = False
        for inds in indices:
            if len(inds) > 0:
                values_masked = True

        if not values_masked:
            print(f"WARNING: No values in rhos match {value}. No values are being masked", file=sys.stderr)

        mask[indices] = mask_value
        return mask

    def mask_water(rhos, water_cond=-1.20397, mask=None, esp=0.00005):
        return ModEMCovariance.mask_value(rhos, water_cond, ModEMCovariance.WATER_MASK, esp=esp)

