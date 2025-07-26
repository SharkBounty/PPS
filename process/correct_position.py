import numpy as np
from read.process import Process
import logging
from process.elevation import Elevation

logging.basicConfig(filename='output.log', level=logging.INFO)

class CorrectPosition():
    def __init__(self, rnx_nav, len_obs, len_sv, tow, sv_health):
        self.sv_health = sv_health
        self.process = Process()
        self.elevation = Elevation()
        self.SPEED_OF_LIGHT = 299792458 
        self.GM = 3.986004418e14
        self.OMEGA_E = 7.2921151467e-5 
        self.F = -2*(np.sqrt(self.GM)/(self.SPEED_OF_LIGHT**2))
        self.rnx_nav = rnx_nav
        self.tow = tow
        self.delta_n = self.process.get_data_avail(self.rnx_nav['DeltaN'], 'G', filter=True, sv_health=self.sv_health)
        self.trans_time = self.process.get_data_avail(self.rnx_nav['TransTime'], 'G', filter=True, sv_health=self.sv_health)
        self.sqrt_a = self.process.get_data_avail(self.rnx_nav['sqrtA'], 'G', filter=True, sv_health=self.sv_health)
        self.m_0 = self.process.get_data_avail(self.rnx_nav['M0'], 'G', filter=True, sv_health=self.sv_health) 
        self.eccentricity = self.process.get_data_avail(self.rnx_nav['Eccentricity'], 'G', filter=True, sv_health=self.sv_health) 
        self.omega = self.process.get_data_avail(self.rnx_nav['omega'], 'G', filter=True, sv_health=self.sv_health)
        self.omega_0 = self.process.get_data_avail(self.rnx_nav['Omega0'], 'G', filter=True, sv_health=self.sv_health) 
        self.omega_dot = self.process.get_data_avail(self.rnx_nav['OmegaDot'], 'G', filter=True, sv_health=self.sv_health) 
        self.cuc = self.process.get_data_avail(self.rnx_nav['Cuc'], 'G', filter=True, sv_health=self.sv_health)
        self.cus = self.process.get_data_avail(self.rnx_nav['Cus'], 'G', filter=True, sv_health=self.sv_health) 
        self.crc = self.process.get_data_avail(self.rnx_nav['Crc'], 'G', filter=True, sv_health=self.sv_health)
        self.crs = self.process.get_data_avail(self.rnx_nav['Crs'], 'G', filter=True, sv_health=self.sv_health) 
        self.cic = self.process.get_data_avail(self.rnx_nav['Cic'], 'G', filter=True, sv_health=self.sv_health) 
        self.cis = self.process.get_data_avail(self.rnx_nav['Cis'], 'G', filter=True, sv_health=self.sv_health) 
        self.io = self.process.get_data_avail(self.rnx_nav['Io'], 'G', filter=True, sv_health=self.sv_health) 
        self.i_dot = self.process.get_data_avail(self.rnx_nav['IDOT'], 'G', filter=True, sv_health=self.sv_health) 
        self.svclockbias = self.process.get_data_avail(self.rnx_nav['SVclockBias'], 'G', filter=True, sv_health=self.sv_health)
        self.svclockdrift = self.process.get_data_avail(self.rnx_nav['SVclockDrift'], 'G', filter=True, sv_health=self.sv_health)
        self.svclockdriftrate = self.process.get_data_avail(self.rnx_nav['SVclockDriftRate'], 'G', filter=True, sv_health=self.sv_health) 
        self.tgd = self.process.get_data_avail(self.rnx_nav['TGD'], 'G', filter=True, sv_health=self.sv_health)
        
        self.dtime = np.array(self.process.get_dtime(self.svclockbias))
        self.toe = self.process.filter(self.rnx_nav.Toe, sv_health)        
        self.corr_ephem = np.full((len_obs, len_sv, 3), np.nan)

        self.Xk = np.full((len_obs, len_sv), np.nan)
        self.Yk = np.full((len_obs, len_sv), np.nan)
        self.Zk = np.full((len_obs, len_sv), np.nan)
        self.dts = np.full((len_obs, len_sv), np.nan)
        self.sat_coord = np.full((len_obs, len_sv, 3), np.nan)
        self.index_tropo = 0
        
    def correct_earth_rotation(self, epoch, obs_time, pseudo_ranges, dtr):
        Xkc = 0
        Ykc = 0
        Zkc = 0

        propagation = (pseudo_ranges/self.SPEED_OF_LIGHT) - dtr + self.dts[epoch]
        
        instant_transmission = self.tow[epoch] - dtr - propagation + self.dts[epoch]
        self.get_correct_pos(obs_time=obs_time, pseudo_ranges=pseudo_ranges, epoch=epoch, instant_transmission=instant_transmission)


        # propagation = propagation[~np.isnan(np.array(propagation))]
        phi = -7.292115e-5 * propagation #7.292115e-5 == self.OMEGA_E (troca só para igualar valores)
        
        for i, iphi in enumerate(phi):
            # Matriz de rotação
                if ~np.isnan(iphi):
                    rotation_matrix = np.array(
                    [
                        [np.cos(iphi), -np.sin(iphi)],
                        [np.sin(iphi), np.cos(iphi)]
                    ]
                    )

                    # Multiplica a matriz de rotação pelas coordenadas efemérides
                    ephem = np.dot(rotation_matrix, [self.Xk[epoch][i], self.Yk[epoch][i]])
                    # self.Xkc.append(ephem[0])
                    # self.Ykc.append(ephem[1])
                    # self.Zkc.append(self.Zk[i])
                    Xkc = ephem[0]
                    Ykc = ephem[1]
                    Zkc = self.Zk[epoch][i]
                else:
                    Xkc = np.nan
                    Ykc = np.nan
                    Zkc = np.nan
                
                self.corr_ephem[epoch][i] = np.column_stack((Xkc, Ykc, Zkc))


        return 

    def ionosphere_correction(self, epoch, X, sat_coord, ion_alpha, ion_beta):
        """
        Calcula o atraso ionosférico baseado no modelo Klobuchar.

        Args:
            tow (float): Tempo da semana (Time of Week).
            X (array): Coordenadas aproximadas do receptor.
            sat_coord (array): Coordenadas do satélite.
            ion_alpha (array): Parâmetros ionosféricos Alpha.
            ion_beta (array): Parâmetros ionosféricos Beta.

        Returns:
            float: Atraso ionosférico em metros.
        """
        rad = np.pi / 180  # Conversão de graus para radianos

        # Calcular coordenadas ENU
        geo = self.elevation.compute_geodesic_coord(*X[0:3])

        # Calcular azimute e elevação
        azimuth, elevation = self.elevation.compute_elevation(X[0:3], sat_coord)

        # Normalizar elevação e azimute
        elevation = rad * elevation / np.pi
        azimuth = rad * azimuth / np.pi

        # Converter latitude e longitude para radianos
        latitude = rad * geo[0] / np.pi
        longitude = rad * geo[1] / np.pi

        # Calcular sub-ponto ionosférico
        psi = 0.0137 / (elevation + 0.11) - 0.022
        sub_iono_lat = np.clip(latitude + psi * np.cos(azimuth), -0.416, 0.416)
        sub_iono_lon = longitude + (psi * np.sin(azimuth) / np.cos(sub_iono_lat * np.pi))

        # Calcular latitude geomagnética
        geomagnetic_lat = sub_iono_lat + 0.064 * np.cos((sub_iono_lon - 1.617) * np.pi)

        # Calcular tempo local no ponto ionosférico
        time_local_ipp = (4.32e4) * sub_iono_lon + self.tow[epoch]
        time_local_ipp = np.where(time_local_ipp > 86400, time_local_ipp - 86400, time_local_ipp)
        time_local_ipp = np.where(time_local_ipp < 0, time_local_ipp + 86400, time_local_ipp)

        # Calcular o fator de inclinação
        slope_factor = 1 + 16 * ((0.53 - elevation) ** 3)

        # Calcular período (P)
        P = ion_beta[0] + ion_beta[1] * geomagnetic_lat + ion_beta[2] * geomagnetic_lat**2 + ion_beta[3] * geomagnetic_lat**3
        P = np.where(P <= 72000, 72000, P)

        # Calcular amplitude (A)
        A = ion_alpha[0] + ion_alpha[1] * geomagnetic_lat + ion_alpha[2] * geomagnetic_lat**2 + ion_alpha[3] * geomagnetic_lat**3
        A = np.where(A < 0, 0, A)

        # Calcular fase (x)
        x = 2 * np.pi * (time_local_ipp - 50400) / P

        # Calcular correção ionosférica
        ic = self.SPEED_OF_LIGHT * slope_factor * (5e-9 + A * (1 - ((x**2) / 2) + ((x**4) / 24)))
        ic = np.where(np.abs(x) > 1.57, self.SPEED_OF_LIGHT * slope_factor * 5e-9, ic)

        return ic


    def troposphere_correction(self, Tc, P, Rh, X, sat_coord):

        rad = np.pi / 180
        azimuth, elevation = self.elevation.compute_elevation(X[0:3], sat_coord)
        
        elevation = rad * elevation
        azimuth = rad * azimuth
        
        T0 = Tc + 273.16
        mh = 1 / np.sin(elevation) + (0.00143 / (np.tan(elevation) + 0.0445))
        mw = 1 / np.sin(elevation) + (0.00035 / (np.tan(elevation) + 0.017))
        
        Hh = 40136 + 148.72 * (T0 - 273.16)
        Hw = 11000
        P0 = 6.1164 * np.power(10, (7.591386 * Tc / (Tc + 240.7263))) * Rh / 100
        
        Tzh = (155.2e-7) * (P / T0) * Hh
        Tzw = (155.2e-7) * (4810 * P0 / (T0 ** 2)) * Hw
        Tzh = np.array(Tzh)
        Tzw = np.array(Tzw)
        
    
        tr = (Tzh * mh) + (Tzw * mw)
        
        return tr

    def get_correct_pos(self, obs_time=None, pseudo_ranges=None, epoch=None, instant_transmission=None):
        print(epoch)
        # self.index = np.argmax(np.where(self.dtime < pseudo_ranges.time.data))
        self.index = np.argmin(abs(obs_time - self.dtime))

        if instant_transmission is None:
            tgps = self.tow[epoch] - (pseudo_ranges/self.SPEED_OF_LIGHT)
            self.pseudo_ranges = pseudo_ranges
        else:
            tgps = instant_transmission
        
        # sv_avail = list(set(self.svclockbias[index].sv.data) & set(pseudo_ranges.sv.data))
        # obs_avail = pseudo_ranges.sel(sv=sv_avail)
        a0 = self.svclockbias[self.index]
        a1 = self.svclockdrift[self.index]
        a2 = self.svclockdriftrate[self.index]
        toe = self.toe[self.index]  # Supondo que `toe` não contém `NaN`
        sqrt_a = self.sqrt_a[self.index]
        delta_n = self.delta_n[self.index]
        m_0 = self.m_0[self.index]
        eccentricity = self.eccentricity[self.index]
        omega = self.omega[self.index]
        cuc = self.cuc[self.index]
        cus = self.cus[self.index]
        crc = self.crc[self.index]
        crs = self.crs[self.index]
        cic = self.cic[self.index]
        cis = self.cis[self.index]
        tgd = self.tgd[self.index]
        io = self.io[self.index]
        i_dot = self.i_dot[self.index]
        omega_0 = self.omega_0[self.index]
        omega_dot = self.omega_dot[self.index]


        n0 = np.sqrt(self.GM/((sqrt_a**2)**3))
        n = n0 + delta_n

        deltaTK = tgps - toe
        deltaTK = np.where(deltaTK > 302400, deltaTK - 604800, deltaTK)
        deltaTK = np.where(deltaTK < -302400, deltaTK + 604800, deltaTK)     

        Mk = m_0 + n * deltaTK

        ek = Mk
        for _ in range(6):
            ek = Mk + eccentricity * np.sin(ek)


        cosVk = (np.cos(ek) - eccentricity) / (1 - (eccentricity * np.cos(ek)))
        senVk = (np.sqrt(1 - eccentricity**2) * np.sin(ek)) / (1 - (eccentricity * np.cos(ek)))

        vk = np.where(((senVk >= 0) & (cosVk >= 0)) | ((senVk < 0) & (cosVk >= 0)), np.arctan(senVk / cosVk), np.arctan(senVk / cosVk) + np.pi)
        
        phik = vk+omega #Argumento de latitude
        deltaUk = (cuc * np.cos(2*phik))+(cus * np.sin(2*phik)) #Correção do argumento de latitutde
        uk = phik + deltaUk #Argumento da latitude corrigido
        phiRk = crc * np.cos(2 * phik)+ crs * np.sin(2*phik) #Correção do raio
        rk = ((sqrt_a**2) * (1-eccentricity * np.cos(ek))) + phiRk #Raio corrigido
        phiIk = (cic*np.cos(2*phik))+(cis*np.sin(2*phik)) #Inclinação corrigida
        ik = io+(i_dot*deltaTK)+phiIk; #Inclinação corrigida

        dtr = self.F * eccentricity * sqrt_a * np.sin(ek)


        #Posição do satélite no plano orbital
        xk = rk*np.cos(uk)
        yk = rk*np.sin(uk)
        
        #Coordenadas terrestres (WGS 84) do satélite
        omegak = omega_0+(omega_dot*deltaTK) - (self.OMEGA_E*tgps)
        self.Xk[epoch] = (xk*np.cos(omegak))-(yk*np.sin(omegak)*np.cos(ik))
        self.Yk[epoch] = (xk*np.sin(omegak))+(yk*np.cos(omegak)*np.cos(ik))
        self.Zk[epoch] = (yk*np.sin(ik))
    
        self.sat_coord[epoch] = np.column_stack((self.Xk[epoch], self.Yk[epoch], self.Zk[epoch]))

        # self.Xk = self.process.clear_nan_values(self.Xk)
        # self.Yk = self.process.clear_nan_values(self.Yk)
        # self.Zk = self.process.clear_nan_values(self.Zk)

        #6109190.065511795
        # print(dtr)
        self.dts[epoch] = a0 + a1 * (tgps - toe) + a2 * (tgps - toe)**2 - tgd + dtr

        


        # return self.Xk, self.Yk, self.Zk, self.dts
