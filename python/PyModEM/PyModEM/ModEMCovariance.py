
from typing import Tuple

from PyModEM import ModEMGrid


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

def write_modem_header(file,):
    file.write(header)
    file.write('\n')

def write_grid_dims(file,
                    nx: int,
                    ny: int,
                    nz: int):
    file.write(f"{nx} {ny} {nz}\n")
    file.write("\n")

def write_smoothing(file, n, smoothing):
    for i in range(0, n, 1):
        file.write(f"{smoothing} ")
    file.write('\n')

def write_modem_smoothing(file,
                          grid: ModEMGrid.ModEMGrid,
                          x_smoothing: float,
                          y_smoothing: float,
                          z_smoothing: float,
                          n_smoothing):
    write_smoothing(file, len(grid.z_dims), x_smoothing)
    write_smoothing(file, len(grid.z_dims), y_smoothing)
    write_smoothing(file, 1, z_smoothing)
    file.write('\n')
    file.write(str(n_smoothing)+'\n')
    file.write('\n')

def write_modem_exception_rules(file, n_rules):
    file.write(str(n_rules)+'\n')

def write_modem_covariance(fname: str,
                          grid: ModEMGrid.ModEMGrid,
                          x_soomthing: float,
                          y_smoothing: float,
                          z_smoothing: float,
                          n_smoothing: int):

    with open(fname, 'w') as file:
        write_modem_header(file)
        write_grid_dims(file, len(grid.x_dims), len(grid.y_dims), len(grid.z_dims))
        write_modem_smoothing(file, grid, x_soomthing, y_smoothing, z_smoothing, n_smoothing)
        write_modem_exception_rules(file, 0)


    

        
        
