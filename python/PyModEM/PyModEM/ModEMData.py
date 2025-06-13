import os
import sys

from typing import Tuple
import numpy as np

from dataclasses import dataclass

NUM_HEADER_LINES = 8

DATA_TYPE_COMPONENT_MAP = {'Full_Impedance' : ['ZXX', 'ZXY', 'ZYX', 'ZYY'],
                         'Off_Diagonal_Impedance' : ['ZXY', 'ZYX'],
                         'Full_Vertical_Components' : ['TX ', 'TY '],
                         'Full_Interstation_TF' : ['MXX', 'MXY', 'MYX', 'MYY'],
                         'Off_Diagonal_Rho_Phase' : ['RHOXY', 'PHSXY', 'RHOYX', 'PHSYX'],
                         'Phase_Tensor' : ['PTXX', 'PTXY', 'PTYX', 'PTYY'],
                         'TE_Impedance' : ['TE'],
                         'TM_Impedance' : ['TM']
                         }

DATA_TYPES_2D = ['TE_Impedance', 'TM_Impedance']
DATA_TYPES_3D = [ 'Full_Impedance',
                 'Off_Diagonal_Impedance',
                 'Full_Vertical_Components',
                 'Full_Interstation_TF',
                 'Off_Diagonal_Rho_Phase',
                 'Phase_Tensor'
                 ]

COMMON_PERIODS = [1.63636e01,
                  2.56e1,
                  5.389474e1,
                  1.024e2,
                  2.155789e2,
                  4.096e2,
                  8.623158e2,
                  1.6384e3,
                  4.681143e3,
                  1.872457e4]

DEFAULT_DATA_DESCRIPTION = "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
DEFAULT_SIGN_CONVENTION = "> exp(+i\omega t)"
DEFAULT_UNITS = "> [V/m]/[T]"


@dataclass
class ModEMDataHeader:
    header: str
    description: str
    data_type: str
    sign_conversion: str
    units: str
    orientation: float
    origin: Tuple[float, float]
    nperiods: int
    nstations: int

@dataclass 
class DataEntry:
    station_code: str
    location_latlon : Tuple[float, float]
    location_xyz: Tuple[float, float, float]
    period : float
    component : str
    real : float
    imag : float
    error : float

@dataclass
class StationData:
    station_code : str
    location_latlon : Tuple[float, float]
    location_xyz : Tuple[float, float, float]
    periods : []

DAT_DATA_FORMAT = "{0:.7f} {1:3} {2:.4f} {3:.3f} {4:.3f} {5:.3f} {6:.3f} {7} {8:.7f} {9:.7f} {10:.8f}"


class Station:
    def __init__(self, station_code=None, latlon=None, xyz=None, data=None):
        self._periods_dict = {}
        if data is not None:
            self.data = data
        elif (station_code is not None and latlon is not None and xyz is not None):
            self.data = StationData(station_code=station_code,
                                    location_latlon=latlon,
                                    location_xyz=xyz,
                                    periods=[])

    def __str__(self):
        return str(f"Station {self.data.station_code} - {len(self.data.periods)} periods")

    def __repr__(self):
        return str(f"<Station {self.data.station_code} - {len(self.data.periods)} periods>")

    @property
    def periods(self):
        return self.data.periods

    @property
    def code(self):
        return self.data.station_code

    def has_period(self, period):
        return period in self.data.periods

    def get_component(self, period : float, component : str):
        for p in self.periods:
            if p.period == period and p.component == component:
                return p

        return None

    def add_period(self, entry: StationData, period):
        self.data.periods.append(period)

    def remove_period(self,  period):
        self.data.periods = [p for p in self.data.periods if (p.data.period != period)]

    def write_data(self, file, components: []):
        for period in self.data.periods:
            if period.data.component in components:
                file.write("{0:.7e} {1:3} {2:.4f} {3:.3f} {4:.3f} {5:.3f} {6:.3f} ".format(
                    period.data.period,
                    self.data.station_code,
                    self.data.location_latlon[0],
                    self.data.location_latlon[1],
                    self.data.location_xyz[0],
                    self.data.location_xyz[1],
                    self.data.location_xyz[2]))
                period.write_data(file)
                file.write("\n")


@dataclass
class PeriodData:
    period : float
    stations : [Station]
    component : str
    real : float
    imag : float
    error : float

class Period:
    def __init__(self, period=None, data=None):
        if data is not None:
            self.data = data
        elif period is not None:
            self.data = PeriodData(period=period,
                                   component=[],
                                   real=None,
                                   imag=None,
                                   error=None,
                                   stations=None)

        
    def __str__(self):
        return str(f"Period:{self.data.period}-{self.data.component}")

    def __repr__(self):
        return str(f"<Period:{self.data.period}-{self.data.component}>")

    @property
    def period(self):
        return self.data.period

    @property
    def component(self):
        return self.data.component
    
    @property
    def real(self):
        return self.data.real

    @property
    def imag(self):
        return self.data.imag

    @property
    def error(self):
        return self.data.error

    @error.setter
    def error(self, value):
        self.data.error = value

    @property
    def stations(self):
        return self.data.stations

    def add_station(self, entry: StationData, station):
        if station not in self.data.stations:
            self.data.stations.append(station)

    def has_station(self, station_code: str) -> bool:
        return station_code in self.data.stations

    def remove_station(self, station_code: str) -> Station:
        return self.data.stations.remove(station_code)

    def write_data(self, file):
        file.write(f"{self.data.component} {self.data.real:.7e} {self.data.imag:.7e} {self.data.error:.7e}") 



class ModEMData():
    def __init__(self, filename=None):
        self._filename = filename

        self._stations = {}
        self._periods = {}
        self.headers = []
        self._components = []
        self._data_types = []
        self._parsing = 0

        if self._filename:
            self.load(self._filename)

    def construct_data_filename(data, comment):
        return f'{comment}.{data.nstations}s.{data.nperiods}p.dat'

    def add_station(self, station: Station):
        self._stations[station.data.station_code] = station 

    def add_period(self, period: float):
        self._periods[period] = period

    @property
    def stations(self):
        return self._stations

    @property
    def nstations(self):
        return len(self._stations)

    def is_new_station(self, station_code: str) -> bool:
        return station_code not in self._stations

    def print_stations(self):
        for _, station in self._stations.items():
            print(station)

    def get_station_latlons(self):
        latlon = []

        for _, station in self._stations.items():
            latlon.append(station.data.location_latlon)

        return np.array(latlon)

    def get_station_xyzs(self):
        xyz = []

        for _, station in self._stations.items():
            xyz.append(station.data.location_xyz)

        return np.array(xyz)

    def remove_station(self, station_code, periods=None):
        station = self._stations[station_code]
        del self._stations[station_code]

        for _, period in self._periods.items():
            if period.has_station(station.code):
                period.remove_station(station.code)

    @property
    def periods(self):
        return self._periods

    @property
    def nperiods(self):
        return len(self._periods)

    def is_new_period(self, period: float) -> bool:
        return period not in self.periods

    def print_periods(self):
        for _, period in self._periods.items():
            print(period)

    def remove_period(self, period, stations=None):
        if period not in self._periods:
            raise ValueError(f"'{period}' is not a valid period for this data")

        period = self._periods[period]
        del self._periods[period.data.period]

        for _, station in self._stations.items():
            station.remove_period(period.data.period)

    @property
    def data_types(self):
        return self._data_types

    @property
    def ndata_types(self):
        return len(self._data_types)

    def add_new_data_type(self, data_type):
        self._data_types += data_type.strip('> ')

    @property
    def component(self):
        return self._component

    @property
    def ncomponents(self):
        return len(self._component)

    def parse_header(self, file) -> ModEMDataHeader:
        line_number = 1

        for line in file:
            if line_number == 1:
                header = line.decode('utf-8').strip()

            if line_number == 2:
                description = line.decode('utf-8').strip()

            if line_number == 3:
                data_type = line.decode('utf-8').strip()

            if line_number == 4:
                sign_conversion = line.decode('utf-8').strip()

            if line_number == 5:
                units = line.decode('utf-8').strip()

            if line_number == 6:
                orientation = float(line.decode('utf-8').strip().split()[1])

            if line_number == 7:
                origin = [float(line.decode('utf-8').strip().split()[1]),
                          float(line.decode('utf-8').strip().split()[2])]

            if line_number == 8:
                nperiods = int(line.decode('utf-8').strip().split()[1])
                nstations = int(line.decode('utf-8').strip().split()[2])

            line_number += 1
            if line_number > NUM_HEADER_LINES:
                break

        return ModEMDataHeader(header = header,
            description = description,
            data_type = data_type,
            sign_conversion = sign_conversion,
            units = units,
            orientation = orientation,
            origin = origin,
            nperiods = nperiods,
            nstations = nstations)

    def get_data_type_component_map(self, header):
        return DATA_TYPE_COMPONENT_MAP[header.data_type.strip('> ')]

    def load(self, filename: str):
        with open(filename, 'rb') as data_file:
            self._parsing = 1

            while (self._parsing > 0):
                header = self.parse_header(data_file)
                self.headers.append(header)

                self.add_new_data_type(header.data_type)

                expected_components = self.get_data_type_component_map(header)

                self._components += expected_components

                n_expected_data =  header.nstations * header.nperiods * len(expected_components)

                self.parse_data(data_file, n_expected_data)
                self._parsing -= 1


    def parse_data(self, file, expected_data: int):
        parsing_data = True
        while(parsing_data):
            line = file.readline().decode('utf-8')

            # We are done finishing
            if line == '': 
                return

            line = line.strip().split()

            # Done reading current data and there is another data type to read in
            if line[0] == '#':
                file.seek(-1, 1) # Seek back one for next header parse
                self._parsing += 1
                return

            period = float(line[0])
            station_station_code = str(line[1])
            gg_lat = float(line[2])
            gg_lon = float(line[3])
            x = float(line[4])
            y = float(line[5])
            z = float(line[6])
            component = str(line[7])
            real = float(line[8])
            imag = float(line[9])
            error = float(line[10])

            entry = DataEntry(
                        station_code=station_station_code,
                        period=period,
                        location_latlon=[gg_lat, gg_lon],
                        location_xyz=[x, y, z],
                        component=component,
                        real=real,
                        imag=imag,
                        error=error
                    )

            self.process_entry(entry)

    def process_entry(self, entry: DataEntry):
        station_code = entry.station_code
        period = entry.period
        component = entry.component

        # If it is a new station create them, or grab the existing one
        if self.is_new_station(station_code):
            station = Station(station_code=entry.station_code,
                              latlon=entry.location_latlon, 
                              xyz=entry.location_xyz)
        else: # Get the existing station
            station = self._stations[station_code]

        period = Period(period=period)
        period.data.component = entry.component
        period.data.real = entry.real
        period.data.imag = entry.imag
        period.data.error = entry.error
        station.add_period(entry, period)

        self._stations[station.data.station_code] = station
        self._periods[period.data.period] = period


    def write_modem_data_header(self, header: ModEMDataHeader, file):
        file.write(f"{header.header}\n")
        file.write(f"{header.description}\n")
        file.write(f"{header.data_type}\n")
        file.write(f"{header.sign_conversion}\n")
        file.write(f"{header.units}\n")
        file.write(f"> {header.orientation}\n")
        file.write(f"> {header.origin[0]} {header.origin[1]}\n")
        file.write(f"> {self.nperiods} {self.nstations}\n")

    def write_data(self, filename: str, comment=None):
        with open(filename, 'w') as data_file:

            for header in self.headers:
                if comment is None:
                    header.header = "# ModEM Data written by PyModEm.ModEMData"
                else:
                    header.header = comment

                # Write the header
                self.write_modem_data_header(header, data_file)
                components_in_header = self.get_data_type_component_map(header)

                for _, station in self.stations.items():
                    station.write_data(data_file, components_in_header)


    def nearest_neighbor(self, station_key, nearest_neighbors=1):
        if station_key not in self.stations:
            raise ValueError(f"'{station_key}' is not in the list of stations")

        station = self.stations[station_key]
        neighbors = [float('inf')] * nearest_neighbors

        for stations in self.stations.values():
            distance = np.linalg.norm(np.array(station.data.location_xyz) - np.array(stations.data.location_xyz))

            if distance == 0.0:
                continue

            index = -1
            for i, dist in enumerate(neighbors):
                if distance < dist:
                    index = i
                    break

            if index != -1:
                neighbors[i] = distance

        return neighbors


def make_period_data_for_datatype(period: float, data_type: str, real=0.0, imag=0.0, error=0.0) -> [Period]:
    if data_type not in DATA_TYPES_3D:
        raise ValueError(f"Data type {data_type} is not a valid data type. Choose from:" \
                          " \'{', '.join(DATA_TYPES_3D}\' ")
    

    components = DATA_TYPE_COMPONENT_MAP[data_type]

    data = []

    for component in components:
        data.append(Period(data=PeriodData(period=period,
                                      stations=[],
                                      component=component,
                                      real=real,
                                      imag=imag,
                                      error=error)))


    return data
    

if __name__ == "__main__":
    import sys
    filename = sys.argv[1]

    data = ModEMData(filename)

    print("Number of stations:", len(data.stations))
    print("Number of periods:", len(data.periods))

    data.write_data("my_test_outfile")
