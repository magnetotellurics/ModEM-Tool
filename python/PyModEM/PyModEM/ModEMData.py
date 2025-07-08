import os
import sys

from typing import List, TextIO, Tuple
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

DEFAULT_DATA_DESCRIPTION = "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
DEFAULT_SIGN_CONVENTION = "> exp(+i\omega t)"
DEFAULT_UNITS = "> [V/m]/[T]"


@dataclass
class ModEMDataHeader:
    ''' ModEMDataHeader - dataclass to hold ModEM data file '''
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
    ''' DataEntry - dataclass to hold a single line in a datafile. Mostly for reading '''
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
    ''' StationData - dataclass to hold station information, including the periods that belong to the that station'''
    station_code : str
    location_latlon : Tuple[float, float]
    location_xyz : Tuple[float, float, float]
    periods : List[any] # A list of PeriodData defined below

''' TODO: Currently we only read/write complex impedeances and full veritical component data types...

need to update code to read/write all MT and CSEM datatypes...
'''
COMPLEX_DATA_FORMAT = "{0:.7f} {1:3} {2:.4f} {3:.3f} {4:.3f} {5:.3f} {6:.3f} {7} {8:.7f} {9:.7f} {10:.8f}"
VERTICAL_DATA_FORMAT = COMPLEX_DATA_FORMAT


class Station:
    ''' Station - Class to handle stations. '''
    def __init__(self, station_code : str = None,
                       latlon       : Tuple[float, float] = None,
                       xyz  : Tuple[float, float, float] = None,
                       data : StationData =None):
        ''' Station(station_code, latlon, xyz, data) - Create a station class 
        
        If `data` is provided, use data as variables for this instatiation otherwise,

        Use station_code, latlon and xyz as new data.

        Parameters
        ----------
        station_code : str
            Station code string name
        latlon : Tuple[float, float]
            Latitude and longitude of the station location
        xyz : Tuple[float, float, float]
            Cartesian coordinates of station location 
        data : StationData
            StationData to use to instatiate this object 
        '''
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

    def has_period(self, period: float) -> bool:
        ''' has_period(period) - Return True if station has data for period,
        otherwise return false.
          
        Parameters
        ----------
        period : float
            Period to query on

        Return
        ------            
        out : bool
            True if station has data for the period, else false.
        '''
        return period in self.data.periods

    def get_component(self, period : float, component : str):
        ''' get_component(period, component) - Get the component for a speicific
        period.
        
        Search the periods associated with this station and return the
        PeriodData that contain the period for the given component. If the
        period and component are not present, return None

        Parameters
        ----------
        period : float
            period to search for
        component : str
            component to search for
        Return
        ------
        out : PeriodData
            The PeriodData that contains the requested period and data, if no periods
        '''
        for p in self.periods:
            if p.period == period and p.component == component:
                return p

        return None

    def add_period(self, period):
        ''' add_period(period) - Add p PeriodData to this station

        Parameters
        ----------
        period : PeriodData
            period to add to this station
        
        '''
        self.data.periods.append(period)

    def remove_period(self,  period):
        ''' remove_period(period) - Remove a period from this station if present

        Parameters
        ----------
        period : PeriodData
            period to remove from this station
        
        '''
        self.data.periods = [p for p in self.data.periods if (p.data.period != period)]

    def write_data(self, file: TextIO, components: List[str]):
        ''' write_data(file, components) - Helper function to write out station data

        This method should be used by the ModEMData class as it just writes out the station data,
        and then calls each of the PeriodDatas to write out their data.

        ModEM Data files can have multiple datatypes in one file. In order to ensure your data 
        is written out correctly, this function should be called with the
        components for only one data type (Although we could update it to just take the datatype
        and use the map above).

        Parameters
        ----------
        file : TextIO
            Open file to write too
        components : List[str]
            List of components to write out. 
        '''
        for period in self.data.periods:
            #TODO: Update this so that instead of components passed in we use the datatype map above
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
    ''' PeriodData - Dataclass to hold a data component associated with a period at a station.

    Note: This dataclass only works for  complex impedances and vertical components. We need to
    update it to use other MT and CSEM datatypes.
    '''
    period : float
    stations : list[Station]
    component : str
    real : float
    imag : float
    error : float

class Period:
    def __init__(self, period : float = None, data : PeriodData = None):
        ''' Period(period, data) - Class to manipulate/hold period data

        If data is present, use data to instansiate this class, if not, create an empty class
        with period set to period.

        Parameters 
        ----------
        period : flaot
            If present (and if data is None), initalize a new class with empty data and period = period
        data : PeriodData
            Use data to instatiate this class
        '''
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

    def add_station(self, station : StationData):
        ''' add_station(station) - Add a station to this period

        Parameters
        ----------
        station : StationData
            Station to add to this period

        '''
        if station not in self.data.stations:
            self.data.stations.append(station)

    def has_station(self, station_code: str) -> bool:
        ''' has_station(station_code) - Return true if period is contains station with code station_code '''
        return station_code in self.data.stations

    def remove_station(self, station_code: str) -> Station:
        ''' remove_station(station_code) - Remove the station associated with station_code from this period '''
        return self.data.stations.remove(station_code)

    def write_data(self, file: TextIO):
        ''' write_data(file) - Write out the data associated with this period 

        This method should be used by the ModEMData class and the Station class
        as it just writes out the station data, and then calls each of the
        PeriodDatas to write out their data.

        Note: This function currently only writes complex impedances and vertical components.
        
        Parameters
        ----------
        file : TextIO
            Open file to write too 
        
        '''
        #TODO: Upadte this function to write non-complex impedances, other MT datatypes and CSEM datatypes
        file.write(f"{self.data.component} {self.data.real:.7e} {self.data.imag:.7e} {self.data.error:.7e}") 



class ModEMData():
    def __init__(self, filename : str = None):
        ''' ModEMData(filename) - Class to read or create a ModEM datafile 

        This class can be used to either read in and manipulate an exisiting ModEM data file, or
        to create a new ModEMData file.

        If filename is present, then read in a data file.
        
        Parameters
        ----------
        filename : str
            If filename is present, read and instatiate this class with the data
            found in filename
        '''
        self._filename = filename

        self._stations = {}
        self._periods = {}
        self.headers = []
        self._components = []
        self._data_types = []
        self._parsing = 0

        if self._filename:
            self.load(self._filename)

    def construct_data_filename(self, comment : str):
        ''' ModEMData.construct_data_filename(comment)

        Create a ModEMData filename in the form of: {comment}.{nstations}.{nperiods}.dat

        Parameter
        ---------
        comment : 
            Comment to append to the front of the filename
        
        '''
        return f'{comment}.{self.nstations}s.{self.nperiods}p.dat'

    def add_station(self, station: Station):
        ''' add_station(station) - Add a new station associated with this Data

        Parameters:
        -----------
        station : Station
            Station to add to this ModEM Data 
        
        '''
        self._stations[station.data.station_code] = station 

    def add_period(self, period: float):
        ''' add_period(period) - Add a new period associated with this Data

        Parameters:
        -----------
        period : float
            Period to add to this ModEM Data 
        
        '''
        self._periods[period] = period

    @property
    def stations(self):
        return self._stations

    @property
    def nstations(self):
        return len(self._stations)

    def is_new_station(self, station_code: str) -> bool:
        ''' is_new_station(station_code) - Return true if station *is not* apart
        of this ModEMData otherwise return false
        '''
        return station_code not in self._stations

    def print_stations(self):
        ''' print_stations - Print all the stations associated with this ModEMData '''
        for _, station in self._stations.items():
            print(station)

    def get_station_latlons(self) -> np.ndarray:
        ''' get_station_latlons() - Return an np.ndarray of all station latitude and longitude '''
        latlon = []

        for _, station in self._stations.items():
            latlon.append(station.data.location_latlon)

        return np.array(latlon)

    def get_station_xyzs(self) -> np.ndarray:
        ''' get_station_xyz() - Return an np.ndarray of all station cartesian locations'''
        xyz = []

        for _, station in self._stations.items():
            xyz.append(station.data.location_xyz)

        return np.array(xyz)

    def remove_station(self, station_code: str):
        ''' remove_station(station_code) - Remove all data entry's for this
        station for this ModEMData'''
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
        ''' is_new_period(period) - Return true if period *is not* apart
        of this ModEMData otherwise return false'''
        return period not in self.periods

    def print_periods(self):
        ''' print_periods() - print all the periods associated with this ModEMData '''
        for _, period in self._periods.items():
            print(period)

    def remove_period(self, period : float):
        ''' remove_period(period) Remove all data entrys that use this period'''
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

    def parse_header(self, file : TextIO) -> ModEMDataHeader:
        ''' parse_header(file) - Read and return a ModEMDataHeader from an open file'''
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

    def get_data_type_component_map(self, header : ModEMDataHeader):
        ''' get_data_type_componet_map(header) -> Get the components associated with this datatype '''
        return DATA_TYPE_COMPONENT_MAP[header.data_type.strip('> ')]

    def load(self, filename: str):
        ''' load(filename : str) - Read and parse a ModEM datafile into this instance'''
        with open(filename, 'rb') as data_file:
            self._parsing = 1

            while (self._parsing > 0):
                # Read the header for this file (or when we encounter a new one)
                header = self.parse_header(data_file)
                self.headers.append(header)

                self.add_new_data_type(header.data_type)
                expected_components = self.get_data_type_component_map(header)
                self._components += expected_components
                n_expected_data =  header.nstations * header.nperiods * len(expected_components)

                # Now, parse the data for the header we read above
                self._parse_data(data_file, n_expected_data)
                self._parsing -= 1


    def _parse_data(self, file : TextIO):
        ''' _parse_data(file) - Read the data associated with a Header '''
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

            # TODO: Update this so that we can read other MT data types and CSEM datatypes

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

            self._process_entry(entry)

    def _process_entry(self, entry: DataEntry):
        ''' _process_entry(entry : DataEntry) -  Process a dentry entry 
        
        This internal function should only be called by _parse_data and that
        should be called by load. This function processes the data entry and
        either adds this data entry to an exisiting station and period, or
        creates a new station and period.
        '''
        station_code = entry.station_code
        period = entry.period

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


    def _write_modem_data_header(self, header: ModEMDataHeader, file):
        ''' _write_modem_data_header(header, file) - Given a ModEMDataHeader object, write out a ModEM header 

        This method should be called from write_data. 
        '''
        file.write(f"{header.header}\n")
        file.write(f"{header.description}\n")
        file.write(f"{header.data_type}\n")
        file.write(f"{header.sign_conversion}\n")
        file.write(f"{header.units}\n")
        file.write(f"> {header.orientation}\n")
        file.write(f"> {header.origin[0]} {header.origin[1]}\n")
        file.write(f"> {self.nperiods} {self.nstations}\n")

    def write_data(self, filename: str, comment=None):
        ''' write_data(filename, comment) - Write out this object in ModEM Data format 

        This function can create a ModEM data file. At this time, you might need to create
        your own header if you create synthetic data. See `make_modem_data` for an example.

        Parameters
        ----------
        filename : str
            Name of file to write too
        comment : str
            Comment to add at the top of the Datafile 
        '''
        with open(filename, 'w') as data_file:

            for header in self.headers:
                if comment is None:
                    header.header = "# ModEM Data written by PyModEm.ModEMData"
                else:
                    header.header = comment

                # Write the header
                self._write_modem_data_header(header, data_file)
                components_in_header = self.get_data_type_component_map(header)

                for _, station in self.stations.items():
                    station.write_data(data_file, components_in_header)


    def nearest_neighbor(self, station_key : str, nearest_neighbors : int = 1):
        ''' nearest_neighbor(station_key, nearest_neighbors=1) - Find the n nearest stations to a station '''
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


def make_period_data_for_datatype(period: float, data_type: str, real=0.0, imag=0.0, error=0.0) -> List[Period]:
    ''' make_period_data_for_datatype(period, data_type, real, imag, error)
    
    Given a period and a data_type, create individual Period class for each component. For instance
    passing in make_period_data_for_datatype(1.0, 'Full_Impedance') returns:

    [<Period:1.0-ZXX>, <Period:1.0-ZXY>, <Period:1.0-ZYX>, <Period:1.0-ZYY>]

    (A list of Periods classes, one for each component)

    Parameters
    ----------
    period : float
        Period to use for the Period class
    data_type : str
        Data type, must match the data type in ModEMData.DATA_TYPE_COMPONENT_MAP
    real : float
        Value to use for the real value, optional 
    imag : float
        value to use for the imaginary part, optional
    error : float
        Error value to use, optional
    
    Returns
    -------
    out : List[Period]
        List of periods, one period for each component.
    '''
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