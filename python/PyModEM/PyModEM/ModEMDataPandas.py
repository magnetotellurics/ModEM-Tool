import io
import sys

import pandas as pd

from PyModEM import ModEMData
from PyModEM.ModEMData import DataEntry

class ModEMDataPandas(ModEMData.ModEMData):

    def __init__(self, fname=None, dataframe=None):
        super().__init__()
        self.dataFrame = None

        if fname is not None:
            self.load(fname)
            self.to_pandas()

        if dataframe is not None and fname is None:
            self.dataFrame = dataframe

    def __sub__(self, other):
        diff_pandas = self.dataFrame - other.dataFrame

        return ModEMDataPandas.from_pandas(diff_pandas, headers=self.headers)

    def from_pandas(df, headers=None):
        data = ModEMDataPandas()

        for index, row in df.iterrows():

            period = index[0]
            station =  index[1]
            component = index[2]

            lat = row['lat']
            lon = row['lon']
            x = row['x']
            y = row['y']
            z = row['z']
            real = row['real']
            imag = row['imag']
            error = row['error']

            entry = DataEntry(
                    station_code = station,
                    location_latlon = (lat, lon),
                    location_xyz = (x, y, z), 
                    period = period,
                    component = component,
                    real = real,
                    imag = imag,
                    error = error
                )

            data._process_entry(entry)

        # Extract the headers for each datatype
        if headers is None:
            for datatype in data.data_types:
                data.headers.append(getattr(df, datatype))
        else:
            data.headers = headers

        data.dataFrame = df

        return data

    def to_pandas(self):
        data = []
        for station in self.stations.values():
            for period in station.periods:
                entries = period.to_entries(station)

                for entry in entries:
                    row = {}

                    row['period'] = entry.period
                    row['station'] = entry.station_code
                    row['component'] = entry.component
                    row['lat'] = entry.location_latlon[0]
                    row['lon'] = entry.location_latlon[1]
                    row['x'] = entry.location_xyz[0]
                    row['y'] = entry.location_xyz[1]
                    row['z'] = entry.location_xyz[2]
                    row['real'] = entry.real
                    row['imag'] = entry.imag
                    row['error'] = entry.error
                    data.append(row)

        self.dataFrame = pd.DataFrame(data)
        self.dataFrame = self.dataFrame.set_index(['period','station','component'])

        # Copy over the headers for later extraction
        for header in self.headers:
            datatype = header.data_type.replace(">","").strip()
            setattr(self.dataFrame, datatype, header)

        return self.dataFrame
