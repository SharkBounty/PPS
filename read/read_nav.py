import numpy as np
import xarray as xr
import pandas as pd
from read.process import Process

class ReadNav():
    def __init__(self, rnx_nav):
        self.rnx_nav = rnx_nav
        self.process = Process()
        self.delta_n = np.array(self.process.get_data_avail(self.rnx_nav['DeltaN'], 'G', filter=True))
        self.trans_time = np.array(self.process.get_data_avail(self.rnx_nav['TransTime'], 'G', filter=True))
        self.sqrt_a = np.array(self.process.get_data_avail(self.rnx_nav['sqrtA'], 'G', filter=True))
        self.m_0 = np.array(self.process.get_data_avail(self.rnx_nav['M0'], 'G', filter=True)) 
        self.eccentricity = np.array(self.process.get_data_avail(self.rnx_nav['Eccentricity'], 'G', filter=True)) 
        self.omega = np.array(self.process.get_data_avail(self.rnx_nav['omega'], 'G', filter=True))
        self.omega_0 = np.array(self.process.get_data_avail(self.rnx_nav['Omega0'], 'G', filter=True)) 
        self.omega_dot = np.array(self.process.get_data_avail(self.rnx_nav['OmegaDot'], 'G', filter=True)) 
        self.cuc = np.array(self.process.get_data_avail(self.rnx_nav['Cuc'], 'G', filter=True))
        self.cus = np.array(self.process.get_data_avail(self.rnx_nav['Cus'], 'G', filter=True)) 
        self.crc = np.array(self.process.get_data_avail(self.rnx_nav['Crc'], 'G', filter=True))
        self.crs = np.array(self.process.get_data_avail(self.rnx_nav['Crs'], 'G', filter=True)) 
        self.cic = np.array(self.process.get_data_avail(self.rnx_nav['Cic'], 'G', filter=True)) 
        self.cis = np.array(self.process.get_data_avail(self.rnx_nav['Cis'], 'G', filter=True)) 
        self.io = np.array(self.process.get_data_avail(self.rnx_nav['Io'], 'G', filter=True)) 
        self.i_dot = np.array(self.process.get_data_avail(self.rnx_nav['IDOT'], 'G', filter=True)) 
        self.svclockbias = self.process.get_data_avail(self.rnx_nav['SVclockBias'], 'G', filter=True)
        self.svclockdrift = self.process.get_data_avail(self.rnx_nav['SVclockDrift'], 'G', filter=True)
        self.svclockdriftrate = self.process.get_data_avail(self.rnx_nav['SVclockDriftRate'], 'G', filter=True) 
        self.dtime = np.array(self.process.get_dtime(self.svclockbias))
        self.toe = self.process.filter(self.rnx_nav.Toe)

    def read(self):
        dataset = xr.Dataset(
            {
                'delta_n': (['time', 'sv'], self.delta_n),
                'trans_time': (['time', 'sv'], self.trans_time),
                'sqrt_a': (['time', 'sv'], self.sqrt_a),
                'm_0': (['time', 'sv'], self.m_0),
                'eccentricity': (['time', 'sv'], self.eccentricity),
                'omega': (['time', 'sv'], self.omega),
                'omega_0': (['time', 'sv'], self.omega_0),
                'omega_dot': (['time', 'sv'], self.omega_dot),
                'cuc': (['time', 'sv'], self.cuc),
                'cus': (['time', 'sv'], self.cus),
                'crc': (['time', 'sv'], self.crc),
                'crs': (['time', 'sv'], self.crs),
                'cic': (['time', 'sv'], self.cic),
                'cis': (['time', 'sv'], self.cis),
                'io': (['time', 'sv'], self.io),
                'i_dot': (['time', 'sv'], self.i_dot),
                'svclockbias': (['time', 'sv'], self.svclockbias),
                'svclockdrift': (['time', 'sv'], self.svclockdrift),
                'svclockdriftrate': (['time', 'sv'], self.svclockdriftrate),
                'toe': (['time', 'sv'], self.toe.data)
            },
            coords={
                'time': self.dtime,
                'sv': self.svclockbias[0].sv.data
            }
        )
        print(dataset)
