import sys
from dataclasses import dataclass
import math
import numpy as np

from ModEMGrid import ModEMGrid, ModEMGridData

def create_idealized_rhos(grid: np.ndarray) -> np.ndarray:
    rhos = grid._data.rhos
    nx = rhos.shape[0]
    ny = rhos.shape[1]
    nz = rhos.shape[2]

    rhos[:,:,:] = 4.60517
    print(nx, ny, nz)

    start_x = nx//3; end_x = nx//2
    start_y = ny//3; end_y = ny - (ny//3)
    print(start_x, end_x, start_y, end_y)
    rhos[start_x:end_x, start_y:end_y, 0:nz//2] = 3.60517

    start_x = nx//2; end_x = (nx - (nx//3))
    start_y = ny//3; end_y = (ny - (ny//3))
    print(start_x, end_x, start_y, end_y)
    rhos[start_x:end_x, start_y:end_y, 0:nz//2] = 6.090776

    return rhos


def make_log_increasing_array(
    z1_layer, target_depth, n_layers, increment_factor=0.9
):
    """Create depth array with log increasing cells, down to target depth,
    inputs are z1_layer thickness, target depth, number of layers (n_layers)
    """

    # make initial guess for maximum cell thickness
    max_cell_thickness = target_depth
    # make initial guess for log_z
    log_z = np.logspace(
        np.log10(z1_layer), np.log10(max_cell_thickness), num=n_layers
    )
    counter = 0

    while np.sum(log_z) > target_depth:
        max_cell_thickness *= increment_factor
        log_z = np.logspace(
            np.log10(z1_layer), np.log10(max_cell_thickness), num=n_layers
        )
        counter += 1
        if counter > 1e6:
            break

    return log_z


def create_mesh_dims(inner_size_e: float,
                inner_size_n: float,
                cell_dim_e: float,
                cell_dim_n: float,
                pad_inner_e : int,
                pad_inner_n : int,
                pad_ext_e : int,
                pad_ext_n : int,
                total_size_e: float,
                total_size_n: float,
                z1_layer=70,
                target_depth=600000,
                nz=60,
                inital_resistivity=1000):

    z = make_log_increasing_array(z1_layer, target_depth, nz)
    print(z.size)

    n_cells_e = int(inner_size_e / cell_dim_e)
    n_cells_n = int(inner_size_n / cell_dim_n)

    n_cells_e += 2 * pad_inner_e
    n_cells_n += 2 * pad_inner_e

    inner_e = np.array([cell_dim_e] * n_cells_e)
    inner_n = np.array([cell_dim_n] * n_cells_n)

    need_total_e = total_size_e - sum(inner_e)
    need_total_n = total_size_n - sum(inner_n)

    need_half_e = need_total_e/2
    need_half_n = need_total_n/2

    extension_e = np.around(make_log_increasing_array(cell_dim_e, need_half_e, pad_ext_e))
    extension_n = np.around(make_log_increasing_array(cell_dim_n, need_half_n, pad_ext_n))

    extension_e_reversed = np.flip(extension_e)
    extension_n_reversed = np.flip(extension_n)

    east_west_dims = np.append(extension_e_reversed, np.append(inner_e, extension_e))
    north_south_dims = np.append(extension_n_reversed, np.append(inner_n, extension_n))

    grid_e = np.append(np.append(-1 * extension_e_reversed + inner_e.min(), inner_e),
                       extension_e + inner_e.max())
    grid_n = np.append(np.append(-1 * extension_n_reversed + inner_n.min(), inner_n),
                       extension_n + inner_n.max())

    center_east = -east_west_dims.__abs__().sum()/1.93
    center_north = -north_south_dims.__abs__().sum()/2.45
    print("Center: ", center_east, center_north)

    rhos = np.ones((east_west_dims.shape[0], north_south_dims.shape[0], z.shape[0]))
    rhos[:,:,:] = math.log(inital_resistivity)
    print(rhos.shape, east_west_dims.shape, north_south_dims.shape, z.shape)

    data = ModEMGridData(
        grid_ew = east_west_dims,
        grid_ns = north_south_dims,
        grid_z = z,
        rhos = rhos,
        origin = [center_east, center_north],
        orientation= 0.0
    )

    return data



inner_size_e = 358400
inner_size_n = 358400
inner_dim_e = 25600
inner_dim_n = 25600
npad_e = 4
npad_n = 4
npad_extension = 10
total_size = 2000000

data = create_mesh_dims(inner_size_e, inner_size_n,
                          inner_dim_e, inner_dim_n,
                          npad_e, npad_n, 
                          npad_extension, npad_extension,
                          total_size, total_size)
grid = ModEMGrid(data=data)
grid.write_grid(f"1000ohms.25.6km.rho")
grid.rhos[:] = np.log(100)
grid.write_grid(f"100ohms.25.6km.rho")

create_idealized_rhos(grid)
grid.write_grid(f"real.25.6km.rho")


inner_dim_e /= 2
inner_dim_n /= 2

data = create_mesh_dims(inner_size_e, inner_size_n,
                          inner_dim_e, inner_dim_n,
                          npad_e, npad_n, 
                          npad_extension, npad_extension,
                          total_size, total_size)
grid = ModEMGrid(data=data)

grid.write_grid(f"1000ohms.12.8km.rho")

grid.rhos[:] = np.log(100)
grid.write_grid(f"100ohms.12.8km.rho")

create_idealized_rhos(grid)
grid.write_grid(f"real.12.8km.rho")

inner_dim_e /= 2
inner_dim_n /= 2

data = create_mesh_dims(inner_size_e, inner_size_n,
                          inner_dim_e, inner_dim_n,
                          npad_e, npad_n, 
                          npad_extension, npad_extension,
                          total_size, total_size)
grid = ModEMGrid(data=data)
grid.write_grid(f"1000ohms.6.4km.rho")
grid.rhos[:] = np.log(100)
grid.write_grid(f"100ohms.6.4km.rho")
create_idealized_rhos(grid)
grid.write_grid(f"real.6.4km.rho")

inner_dim_e /= 2
inner_dim_n /= 2

data = create_mesh_dims(inner_size_e, inner_size_n,
                          inner_dim_e, inner_dim_n,
                          npad_e, npad_n, 
                          npad_extension, npad_extension,
                          total_size, total_size)
grid = ModEMGrid(data=data)
grid.write_grid(f"1000ohms.3.2km.rho")
grid.rhos[:] = np.log(100)
grid.write_grid(f"100ohms.3.2km.rho")
create_idealized_rhos(grid)
grid.write_grid(f"real.32km.rho")

inner_dim_e /= 2
inner_dim_n /= 2

data = create_mesh_dims(inner_size_e, inner_size_n,
                          inner_dim_e, inner_dim_n,
                          npad_e, npad_n, 
                          npad_extension, npad_extension,
                          total_size, total_size)
grid = ModEMGrid(data=data)
grid.write_grid(f"1000ohms.1.6km.rho")
grid.rhos[:] = np.log(100)
grid.write_grid(f"100ohms.1.6km.rho")
create_idealized_rhos(grid)
grid.write_grid(f"real.1.6km.rho")

inner_dim_e /= 2
inner_dim_n /= 2

data = create_mesh_dims(inner_size_e, inner_size_n,
                          inner_dim_e, inner_dim_n,
                          npad_e, npad_n, 
                          npad_extension, npad_extension,
                          total_size, total_size)
grid = ModEMGrid(data=data)
grid.write_grid(f"1000ohms..8km.rho")
grid.rhos[:] = np.log(100)
grid.write_grid(f"100ohms..8km.rho")
create_idealized_rhos(grid)
grid.write_grid(f"real..8km.rho")
