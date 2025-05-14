import os

from typing import Tuple

import numpy as np

def parse_mackie_header(header: str) -> [int, int, str]:
    header_split = header.strip().split()

    if len(header_split) != 3:
        raise ValueError(f"{filename} does not appear to be a 2D Resistivity model in Randi Mackie's format.")

    try:
        ny = int(header_split[0])
    except ValueError as error:
        raise ValueError(f"ny {header_split[0]} needs to be an integer value") from error

    try:
        nz = int(header_split[1])
    except ValueError as error:
        raise ValueError(f"nz {header_split[1]} needs to be an integer value") from error

    resistivity_type = header_split[2]

    if resistivity_type != 'LOGE' and resistivity_type != 'LINEAR':
        raise ValueError(f"Error: Unepexcted resistivity type: {resistivity_type}")

    return ny, nz, resistivity_type 


def read_header(filename: str) -> Tuple[int, int, str]:
    with open(filename, 'r') as model_file:
        header = model_file.readline()
        ny, nz, resistivity_type = parse_mackie_header(header)

    return ny, nz, resistivity_type 

def read_dim(file, n):
    count = 0
    arry = []
    while count != n:
        line = file.readline().strip().split()
        count += len(line)
        arry += line

    arry = np.array(arry)
    return arry.astype(float)

def read_rhos(file, ny, nz):
    count = 0
    rhos = []

    while count != ny * nz:
        line = file.readline().strip().split()
        count += len(line)
        rhos += line

    rhos = np.array(rhos)
    rhos = rhos.astype(float)
    return np.reshape(rhos, (nz, ny))

def load2(filename: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    with open(filename, 'r') as model_file:
        header = model_file.readline()

        ny, nz, resistivity_type = parse_mackie_header(header)

        y = read_dim(model_file, ny)
        z = read_dim(model_file, nz)
        
        tmp = model_file.readline()

        rho = read_rhos(model_file, ny, nz)


        return y, z, rho, resistivity_type


def load(filename: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    with open(filename, 'r') as model_file:
        header = model_file.readline()

        ny, nz, resistivity_type = parse_mackie_header(header)

        y = np.fromfile(model_file, count=ny)
        z = np.fromfile(model_file, dtype='f', count=nz)
        tmp = np.fromfile(model_file, dtype='f', count=1)
        rhos = np.fromfile(model_file, dtype='f', count=ny * nz)

        rhos = np.reshape(rhos, (ny, nz))

        return y, z, rhos, resistivity_type
        
if __name__ == "__main__":
    import sys
    filename = sys.argv[1]
    y, z, rhos, res_type = load(filename)
    print(y.shape, z.shape, rhos.shape, res_type)
    print(read_header(filename))
    print(rhos)
