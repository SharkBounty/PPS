# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

import numpy as np
cimport numpy as np
from libc.math cimport sin, cos, sqrt, M_PI, isnan, tan, floor, abs, atan2


# Importa as classes Python originais e xarray
from read.process import Process
from process.elevation import Elevation
import xarray as xr

# cdef class define a classe como um tipo de extensão C para performance máxima
cdef class CorrectPosition:
    # --- DECLARAÇÃO DE ATRIBUTOS ---
    # Atributos de dados de navegação são 'object' para guardar objetos xarray com metadados
    cdef public object sv_health, process, elevation, rnx_nav, tow, pseudo_ranges
    cdef public object delta_n, trans_time, sqrt_a, m_0, eccentricity, omega, omega_0, omega_dot
    cdef public object cuc, cus, crc, crs, cic, cis, io, i_dot
    cdef public object svclockbias, svclockdrift, svclockdriftrate, tgd, toe
    cdef public object dtime # Mantido como object para flexibilidade

    # Atributos de resultado são 'np.ndarray' para performance nos cálculos
    cdef public np.ndarray corr_ephem, Xk, Yk, Zk, dts, sat_coord
    cdef public int index, index_tropo
    
    # Constantes
    cdef public double SPEED_OF_LIGHT, GM, OMEGA_E, F

    # --- MÉTODO __INIT__ ---
    # Projetado para trabalhar com o objeto rnx_nav já unificado
    def __init__(self, rnx_nav, len_obs, len_sv, tow, sv_health):
        self.sv_health = sv_health
        self.process = Process()
        self.elevation = Elevation()

        # Constantes
        self.SPEED_OF_LIGHT = 299792458.0
        self.GM = 3.986004418e14
        self.OMEGA_E = 7.2921151467e-5
        self.F = -2 * (sqrt(self.GM) / (self.SPEED_OF_LIGHT**2))

        # Trabalha diretamente com os objetos recebidos
        self.rnx_nav = rnx_nav
        self.tow = tow

        # Extrai os dados do objeto rnx_nav principal
        self.delta_n = self.process.get_data_avail(rnx_nav['DeltaN'], 'G', filter=True, sv_health=sv_health)
        self.trans_time = self.process.get_data_avail(rnx_nav['TransTime'], 'G', filter=True, sv_health=sv_health)
        self.sqrt_a = self.process.get_data_avail(rnx_nav['sqrtA'], 'G', filter=True, sv_health=sv_health)
        self.m_0 = self.process.get_data_avail(rnx_nav['M0'], 'G', filter=True, sv_health=sv_health)
        self.eccentricity = self.process.get_data_avail(rnx_nav['Eccentricity'], 'G', filter=True, sv_health=sv_health)
        self.omega = self.process.get_data_avail(rnx_nav['omega'], 'G', filter=True, sv_health=sv_health)
        self.omega_0 = self.process.get_data_avail(rnx_nav['Omega0'], 'G', filter=True, sv_health=sv_health)
        self.omega_dot = self.process.get_data_avail(rnx_nav['OmegaDot'], 'G', filter=True, sv_health=sv_health)
        self.cuc = self.process.get_data_avail(rnx_nav['Cuc'], 'G', filter=True, sv_health=sv_health)
        self.cus = self.process.get_data_avail(rnx_nav['Cus'], 'G', filter=True, sv_health=sv_health)
        self.crc = self.process.get_data_avail(rnx_nav['Crc'], 'G', filter=True, sv_health=sv_health)
        self.crs = self.process.get_data_avail(rnx_nav['Crs'], 'G', filter=True, sv_health=sv_health)
        self.cic = self.process.get_data_avail(rnx_nav['Cic'], 'G', filter=True, sv_health=sv_health)
        self.cis = self.process.get_data_avail(rnx_nav['Cis'], 'G', filter=True, sv_health=sv_health)
        self.io = self.process.get_data_avail(rnx_nav['Io'], 'G', filter=True, sv_health=sv_health)
        self.i_dot = self.process.get_data_avail(rnx_nav['IDOT'], 'G', filter=True, sv_health=sv_health)
        self.svclockbias = self.process.get_data_avail(rnx_nav['SVclockBias'], 'G', filter=True, sv_health=sv_health)
        self.svclockdrift = self.process.get_data_avail(rnx_nav['SVclockDrift'], 'G', filter=True, sv_health=sv_health)
        self.svclockdriftrate = self.process.get_data_avail(rnx_nav['SVclockDriftRate'], 'G', filter=True, sv_health=sv_health)
        self.tgd = self.process.get_data_avail(rnx_nav['TGD'], 'G', filter=True, sv_health=sv_health)
        self.toe = self.process.filter(rnx_nav.Toe, sv_health)

        # A função get_dtime recebe o objeto xarray completo e deve funcionar
        self.dtime = np.array(self.process.get_dtime(self.svclockbias))
        
        # Inicialização dos arrays de resultados
        self.corr_ephem = np.full((len_obs, len_sv, 3), np.nan, dtype=np.double)
        self.Xk = np.full((len_obs, len_sv), np.nan, dtype=np.double)
        self.Yk = np.full((len_obs, len_sv), np.nan, dtype=np.double)
        self.Zk = np.full((len_obs, len_sv), np.nan, dtype=np.double)
        self.dts = np.full((len_obs, len_sv), np.nan, dtype=np.double)
        self.sat_coord = np.full((len_obs, len_sv, 3), np.nan, dtype=np.double)
        self.index_tropo = 0
        self.index = 0

    # --- MÉTODOS DE CÁLCULO OTIMIZADOS ---
    cpdef get_correct_pos(self, obs_time=None, pseudo_ranges=None, epoch=None, instant_transmission=None):
        # Linha corrigida
        self.index = np.argmin(np.abs(obs_time - self.dtime))
        
        
        cdef np.ndarray[np.double_t, ndim=1] tgps_arr
        if instant_transmission is None:
            tgps_arr = self.tow[epoch] - (pseudo_ranges.data / self.SPEED_OF_LIGHT)
            self.pseudo_ranges = pseudo_ranges
        else:
            tgps_arr = instant_transmission
        
        # O acesso com .data extrai o array NumPy do objeto xarray
        cdef np.ndarray[np.double_t, ndim=1] a0 = self.svclockbias[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] a1 = self.svclockdrift[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] a2 = self.svclockdriftrate[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] toe = self.toe[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] sqrt_a = self.sqrt_a[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] delta_n = self.delta_n[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] m_0 = self.m_0[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] eccentricity = self.eccentricity[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] omega = self.omega[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] cuc = self.cuc[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] cus = self.cus[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] crc = self.crc[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] crs = self.crs[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] cic = self.cic[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] cis = self.cis[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] tgd = self.tgd[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] io = self.io[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] i_dot = self.i_dot[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] omega_0 = self.omega_0[self.index].data
        cdef np.ndarray[np.double_t, ndim=1] omega_dot = self.omega_dot[self.index].data

        cdef int n_sats = tgps_arr.shape[0]
        cdef np.ndarray[np.double_t, ndim=1] Xk_epoch = self.Xk[epoch]
        cdef np.ndarray[np.double_t, ndim=1] Yk_epoch = self.Yk[epoch]
        cdef np.ndarray[np.double_t, ndim=1] Zk_epoch = self.Zk[epoch]
        cdef np.ndarray[np.double_t, ndim=1] dts_epoch = self.dts[epoch]

        cdef int i, j
        cdef double n0, n, deltaTK, Mk, ek, cosVk, senVk, vk
        cdef double phik, deltaUk, uk, phiRk, rk, phiIk, ik, dtr, xk, yk, omegak

        for i in range(n_sats):
            if isnan(tgps_arr[i]): continue
            
            n0 = sqrt(self.GM / ((sqrt_a[i]**2)**3))
            n = n0 + delta_n[i]
            deltaTK = tgps_arr[i] - toe[i]
            if deltaTK > 302400.0: deltaTK -= 604800.0
            if deltaTK < -302400.0: deltaTK += 604800.0
            Mk = m_0[i] + n * deltaTK
            ek = Mk
            for j in range(6):
                ek = Mk + eccentricity[i] * sin(ek)
            
            cosVk = (cos(ek) - eccentricity[i]) / (1.0 - (eccentricity[i] * cos(ek)))
            senVk = (sqrt(1.0 - eccentricity[i]**2) * sin(ek)) / (1.0 - (eccentricity[i] * cos(ek)))
            
            vk = atan2(senVk, cosVk)
            
            phik = vk + omega[i]
            deltaUk = (cuc[i] * cos(2 * phik)) + (cus[i] * sin(2 * phik))
            uk = phik + deltaUk
            phiRk = crc[i] * cos(2 * phik) + crs[i] * sin(2 * phik)
            rk = ((sqrt_a[i]**2) * (1.0 - eccentricity[i] * cos(ek))) + phiRk
            phiIk = (cic[i] * cos(2 * phik)) + (cis[i] * sin(2 * phik))
            ik = io[i] + (i_dot[i] * deltaTK) + phiIk
            dtr = self.F * eccentricity[i] * sqrt_a[i] * sin(ek)
            xk = rk * cos(uk)
            yk = rk * sin(uk)
            omegak = omega_0[i] + (omega_dot[i] * deltaTK) - (self.OMEGA_E * tgps_arr[i])
            
            Xk_epoch[i] = (xk * cos(omegak)) - (yk * sin(omegak) * cos(ik))
            Yk_epoch[i] = (xk * sin(omegak)) + (yk * cos(omegak) * cos(ik))
            Zk_epoch[i] = yk * sin(ik)
            dts_epoch[i] = a0[i] + a1[i] * deltaTK + a2[i] * deltaTK**2 - tgd[i] + dtr

        self.sat_coord[epoch] = np.column_stack((Xk_epoch, Yk_epoch, Zk_epoch))

    cpdef correct_earth_rotation(self, int epoch, object obs_time, np.ndarray[np.double_t, ndim=1] pseudo_ranges, double dtr_clock_error):
        # Acessa os dados .data para garantir que são arrays NumPy
        cdef np.ndarray[np.double_t, ndim=1] dts_epoch = self.dts[epoch]
        cdef np.ndarray[np.double_t, ndim=1] propagation = (pseudo_ranges / self.SPEED_OF_LIGHT) - dtr_clock_error + dts_epoch
        cdef np.ndarray[np.double_t, ndim=1] instant_transmission = self.tow[epoch] - dtr_clock_error - propagation + dts_epoch
        
        self.get_correct_pos(obs_time=obs_time, pseudo_ranges=pseudo_ranges, epoch=epoch, instant_transmission=instant_transmission)

        cdef np.ndarray[np.double_t, ndim=1] phi = -self.OMEGA_E * propagation
        cdef int n_sats = phi.shape[0]
        cdef int i
        cdef double iphi, c_phi, s_phi, ephem0, ephem1
        
        cdef np.ndarray[np.double_t, ndim=1] Xk_epoch = self.Xk[epoch]
        cdef np.ndarray[np.double_t, ndim=1] Yk_epoch = self.Yk[epoch]
        cdef np.ndarray[np.double_t, ndim=1] Zk_epoch = self.Zk[epoch]

        for i in range(n_sats):
            iphi = phi[i]
            if not isnan(iphi):
                c_phi = cos(iphi)
                s_phi = sin(iphi)
                ephem0 = c_phi * Xk_epoch[i] - s_phi * Yk_epoch[i]
                ephem1 = s_phi * Xk_epoch[i] + c_phi * Yk_epoch[i]
                
                self.corr_ephem[epoch, i, 0] = ephem0
                self.corr_ephem[epoch, i, 1] = ephem1
                self.corr_ephem[epoch, i, 2] = Zk_epoch[i]
            else:
                self.corr_ephem[epoch, i, 0] = np.nan
                self.corr_ephem[epoch, i, 1] = np.nan
                self.corr_ephem[epoch, i, 2] = np.nan
        return

    cpdef ionosphere_correction(self, int epoch, np.ndarray[np.double_t, ndim=1] X, np.ndarray[np.double_t, ndim=2] sat_coord, 
                               np.ndarray[np.double_t, ndim=1] ion_alpha, np.ndarray[np.double_t, ndim=1] ion_beta):

        geo = self.elevation.compute_geodesic_coord(X[0], X[1], X[2])
        azimuth_py, elevation_py = self.elevation.compute_elevation(X[0:3], sat_coord)
        
        cdef np.ndarray[np.double_t, ndim=1] azimuth = np.asarray(azimuth_py, dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] elevation = np.asarray(elevation_py, dtype=np.double)

        cdef int n_sats = azimuth.shape[0]
        cdef np.ndarray[np.double_t, ndim=1] ic = np.empty(n_sats, dtype=np.double)
        cdef double lat_rad = geo[0] * M_PI / 180.0
        cdef double lon_rad = geo[1] * M_PI / 180.0
        cdef double current_tow = self.tow[epoch]
        
        cdef int i
        cdef double elev_rad, azim_rad, psi, sub_iono_lat, sub_iono_lon, geo_lat, time_local, sf, P, A, x_phase

        for i in range(n_sats):
            if isnan(elevation[i]):
                ic[i] = np.nan
                continue

            elev_rad = elevation[i] * M_PI / 180.0
            azim_rad = azimuth[i] * M_PI / 180.0
            
            psi = 0.0137 / (elev_rad/M_PI + 0.11) - 0.022
            sub_iono_lat = lat_rad/M_PI + psi * cos(azim_rad)
            if sub_iono_lat > 0.416: sub_iono_lat = 0.416
            if sub_iono_lat < -0.416: sub_iono_lat = -0.416
            
            sub_iono_lon = lon_rad/M_PI + (psi * sin(azim_rad) / cos(sub_iono_lat * M_PI))
            geo_lat = sub_iono_lat + 0.064 * cos((sub_iono_lon - 1.617) * M_PI)
            
            time_local = 4.32e4 * sub_iono_lon + current_tow
            time_local = time_local - floor(time_local / 86400.0) * 86400.0
            if time_local < 0: time_local += 86400.0

            sf = 1.0 + 16.0 * (0.53 - elev_rad/M_PI)**3
            
            P = ion_beta[0] + ion_beta[1] * geo_lat + ion_beta[2] * geo_lat**2 + ion_beta[3] * geo_lat**3
            if P < 72000.0: P = 72000.0

            A = ion_alpha[0] + ion_alpha[1] * geo_lat + ion_alpha[2] * geo_lat**2 + ion_alpha[3] * geo_lat**3
            if A < 0.0: A = 0.0

            x_phase = 2 * M_PI * (time_local - 50400.0) / P

            if abs(x_phase) > 1.57:
                ic[i] = self.SPEED_OF_LIGHT * sf * 5e-9
            else:
                ic[i] = self.SPEED_OF_LIGHT * sf * (5e-9 + A * (1.0 - (x_phase**2 / 2.0) + (x_phase**4 / 24.0)))
        return ic

    cpdef troposphere_correction(self, double Tc, double P, double Rh, np.ndarray[np.double_t, ndim=1] X, np.ndarray[np.double_t, ndim=2] sat_coord):
        azimuth_py, elevation_py = self.elevation.compute_elevation(X[0:3], sat_coord)
        cdef np.ndarray[np.double_t, ndim=1] elevation = np.asarray(elevation_py, dtype=np.double)
        
        cdef int n_sats = elevation.shape[0]
        cdef np.ndarray[np.double_t, ndim=1] tr = np.empty(n_sats, dtype=np.double)
        cdef int i
        cdef double elev_rad, T0, mh, mw, Hh, Hw, P0, Tzh, Tzw
        
        T0 = Tc + 273.16
        Hh = 40136.0 + 148.72 * (T0 - 273.16)
        Hw = 11000.0
        P0 = 6.1164 * (10**(7.591386 * Tc / (Tc + 240.7263))) * Rh / 100.0
        Tzh = (155.2e-7) * (P / T0) * Hh
        Tzw = (155.2e-7) * (4810.0 * P0 / (T0**2)) * Hw

        for i in range(n_sats):
            if isnan(elevation[i]):
                tr[i] = np.nan
                continue

            elev_rad = elevation[i] * M_PI / 180.0
            mh = 1.0 / sin(elev_rad) + (0.00143 / (tan(elev_rad) + 0.0445))
            mw = 1.0 / sin(elev_rad) + (0.00035 / (tan(elev_rad) + 0.017))
            
            tr[i] = (Tzh * mh) + (Tzw * mw)
            
        return tr