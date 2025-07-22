import numpy as np
from read.process import Process
import xarray as xr
from xarray.groupers import UniqueGrouper as ug
import warnings
import georinex as rnx
import os
import time
import pyinstrument
from process.correct_position import CorrectPosition
from process.elevation import Elevation
from scipy.stats import chi2
from matplotlib import pyplot as plt 
import pickle
import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning, message="Converting non-nanosecond.*")

class SPP():
    #Definição das constantes
     

    def __init__(self, paths, start, end, type_obs):
        self.process = Process()
        self.convert_enu = Elevation()

        self.deg_elevation_mask = 0
        self.epoch = 0 
        self.tgps = []
        self.path = paths
        self.rnx_nav = []
        self.rnx_obs = []
        self.start = self.process.date_to_doy(*start)
        self.end = self.process.date_to_doy(*end)
        self.tow = self.process.get_tow(*start, *end)
        self.__read_files()
        self.ionBeta = self.rnx_nav.ionospheric_corr_GPS[4:8]
        self.ionAlpha = self.rnx_nav.ionospheric_corr_GPS[0:4]
        
        data = self.rnx_obs[type_obs]

        self.error_3d = np.zeros(len(data))


        # self.rnx_nav = Rinex(path_files['data'])  
        # self.rnx_obs.load_data(path_files['obs'])
        # self.rnx_nav.load_data(path_files['nav'])
        data = self.rnx_obs[type_obs]
        self.obs_time = self.rnx_obs['time']
        
        self.approx_position = self.rnx_obs.position
        self.brdc_approx_pos = np.array([self.approx_position[0], self.approx_position[1], self.approx_position[2]])

        self.pos_calc = np.zeros((len(data), 3))
        

        self.avail_obs = self.process.get_data_avail(data, 'G')
        self.len_sv = len(self.avail_obs[0].sv)
        self.cpos = CorrectPosition(self.rnx_nav, len(self.avail_obs), self.len_sv, self.tow, self.rnx_nav.health)
        
        self.elevation = np.full( ((len(self.avail_obs), self.len_sv)), np.nan)
        self.pseudo_ranges = np.full(((len(self.avail_obs), self.len_sv)), np.nan)
        # self.mask = np.full(((len(self.avail_obs), self.len_sv)), np.nan)
        self.mask = []

    def update_position(self, A, P, dX, dL, c, Lb, sigLb, sigma0_threshold=0.05, value_test=2.70554345409542, max_repet=10):
        repet = 0
        test = 0
        
        # Passo de atualização da posição
        Qx = np.linalg.inv(A.T @ P @ A)  # Matriz de covariância dos parâmetros ajustados
        dX[3] = dX[3] * c  # Ajusta a quarta componente de dX pela velocidade da luz
        V = A @ dX + dL  # Calcula os resíduos ajustados

        # Calcular sigma0pos
        sigma0pos = (V.T @ P @ V) / (len(Lb) - 4)
        alpha = sigma0_threshold  # Nível de significância para o teste estatístico
        
        # Teste estatístico usando a distribuição qui-quadrado
        if sigma0pos > chi2.ppf(1 - alpha, len(Lb) - 4):
            repet += 1
            test = 1
        else:
            repet = 0
            test = 0

        CovX = sigma0pos * Qx  # Matriz de covariância escalonada
        sigLa = A @ CovX @ A.T  # Matriz de covariância das observações ajustadas
        Qv = sigLb + sigLa  # Matriz de covariância das observações
        Vpad = V / np.sqrt(np.diag(Qv))  # Resíduos padronizados

        # Ajustes baseados em testes estatísticos
        if test == 1:
            if np.max(np.abs(Vpad)) > value_test:
                condition = np.where(np.abs(Vpad) > value_test)[0]
                # Aqui você pode atualizar os dados de entrada conforme necessário
                # Exemplo:
                # posInputs['epoch'][iEpoch]['GPSSat']['Coord'][condition, 4] = 0
            else:
                repet = 0

        if repet == max_repet:
            repet = 0
        return repet
    # return V, CovX, sigma0pos, Qx, sigLa, Qv, Vpad, test, repet


 
    def __read_met(self, file_path):
        # Ignorar as linhas de cabeçalho até o final do cabeçalho
        print(file_path)
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                if 'END OF HEADER' in line:
                    header_line = i + 1
                    break

        # Ler o arquivo como um DataFrame, ignorando as linhas de cabeçalho
        df = pd.read_csv(file_path, delim_whitespace=True, skiprows=header_line, header=None)

        # Renomear as colunas com base nas informações dos dados
        df.columns = ['Ano', 'Mes', 'Dia', 'Hora', 'Minuto', 'Segundo', 'Pressao', 'Temperatura', 'Umidade', 'Velocidade_Vento', 'Direcao_Vento', 'Chuva', 'Insolacao']
        da = xr.DataArray(df.values, dims=('time', 'vars'))
        return da



    def __read_files(self):
        nav = []
        obs = []
        met = []
        for dir in self.path:
            for filename in os.listdir(dir):
                if(filename.endswith('n')):
                    nav.append(rnx.load(os.path.join(dir, filename)))
                elif filename.endswith('o') and self.start <= int(dir.split('/')[-1]) <= self.end:
                    obs.append(rnx.load(os.path.join(dir, filename)))
                elif filename.endswith('m'): 
                    met.append(self.__read_met(os.path.join(dir, filename)))
        self.rnx_nav = xr.concat(nav, dim='time')
        self.rnx_obs = xr.concat(obs, dim='time')
        self.rnx_met = xr.concat(met, dim='time')
            

    
    def __update(self, pseudo_ranges, epoch):
        if epoch == 839:
            print('a')

        mask = self.mask[epoch]


        self.x = np.array([
            self.approx_position[0],
            self.approx_position[1],
            self.approx_position[2],
            0
        ])
                
        elevation = self.elevation[epoch]
        elevation = self.process.clear_nan_values(elevation)
        

        # while self.repet >= 1:
        elevation = elevation/elevation.max()

        P = np.zeros((elevation.size, elevation.size))

        # Inicializa sigLb como uma matriz de zeros do mesmo tamanho
        sigLb = np.zeros_like(P)
        sig0 = 2

        # Loop para calcular sigLb e P
        for i in range(elevation.size):
            sigLb[i, i] = sig0**2 / elevation[i]**2
            P[i, i] = 1 / sigLb[i, i]

        Lb = pseudo_ranges
        error = 1
        j = 1
        index = 0
        while error > 0.0001:
            if (epoch % 4 == 0) and epoch !=0:
                index += 1
            Pr  = self.rnx_met[index][6]
            Tc = self.rnx_met[index][7]
            Rh = self.rnx_met[index][8]

            self.cpos.correct_earth_rotation(epoch, self.obs_time[epoch].data, pseudo_ranges, self.x[3])
            corr_ephem = self.cpos.corr_ephem[epoch]
            corr_ephem[mask] = np.nan
            self.cpos.dts[mask] = np.nan

            T = self.cpos.troposphere_correction(Tc, Pr, Rh, self.x, corr_ephem)
            I = self.cpos.ionosphere_correction(epoch, self.x, corr_ephem, self.ionAlpha, self.ionBeta)

            xaux = (self.x[0] - corr_ephem[:, 0])**2
            yaux = (self.x[1] - corr_ephem[:, 1])**2
            zaux = (self.x[2] - corr_ephem[:, 2])**2

            geom_dist = np.sqrt(xaux + yaux + zaux)
    
            Lo = geom_dist + self.cpos.SPEED_OF_LIGHT * ((self.x[3] - self.cpos.dts[epoch]))
            dL = np.array(Lb) - np.array(Lo)
            dL = self.process.clear_nan_values(dL)

            corr_ephem_filtered = self.process.clear_nan_values(corr_ephem)
            geom_dist_filtered = self.process.clear_nan_values(geom_dist)

            A = np.ones((len(geom_dist_filtered), 4))
            for i in range(3): 
                A[:, i] = (self.x[i] - corr_ephem_filtered[:, i]) / geom_dist_filtered
            
            try:
                inv_ATA = np.linalg.inv((A.T @ P) @ A)
            except np.linalg.LinAlgError:
                inv_ATA = np.linalg.pinv((A.T @ P) @ A)
            dX = inv_ATA @ (A.T @ P) @ dL 
            dX[3] = dX[3] / self.cpos.SPEED_OF_LIGHT
            self.x += dX 
            error = np.linalg.norm(dX)
            

            
            if j > 10:
                print('foi')
                self.x = np.ones(3)*np.nan
                break
            
            j += 1 
            

            # self.repet = self.update_position(A, P, dX, dL, self.cpos.SPEED_OF_LIGHT, Lb, sigLb)    
            # print(repet)
        
        
        self.pos_calc[epoch] = self.x[0:3]
        # print(np.sqrt((self.approx_position[0]-self.x[0])**2+(self.approx_position[1]-self.x[1])**2+(self.approx_position[2]-self.x[2])**2))
        # self.error_3d[epoch] = error_3d   
            
            # print(np.linalg.norm(self.approx_position))
        #2.037457783927498e+07
        # sv_avail = list(set(self.svclockbias[index].sv.data) & set(pseudo_ranges.sv.data))
        # obs_avail = pseudo_ranges.sel(sv=sv_avail)


            

# with open("exemplo.txt", "w") as file_txt:
#     file_txt.write("Olá, mundo!\n")
#     file_txt.write("Este é um exemplo de escrita em um arquivo de texto.\n")

    def run(self):
        self.repet = 1
        # print(self.avail_obs[0])
        # with pyinstrument.profile():
        start = time.perf_counter()

        for epoch, obs in enumerate(self.avail_obs):
            sv = obs.sv
            self.cpos.get_correct_pos(obs_time=obs.time.data, pseudo_ranges=obs, epoch=epoch)
            
            _, self.elevation[epoch] = self.convert_enu.compute_elevation(self.brdc_approx_pos, self.cpos.sat_coord[epoch])

            self.mask.append(np.where(self.elevation[epoch] <= self.deg_elevation_mask)[0])
            
            if len(np.where(self.elevation[epoch] <= self.deg_elevation_mask)[0]) > 0:

                mask = self.mask[epoch]
                obs[mask] = np.nan
                self.elevation[epoch][mask] = np.nan

            self.pseudo_ranges[epoch] = obs

            # if len(mask_indexs) > 0:
            #     self.mask[epoch, 0:len(mask_indexs)] = mask_indexs
            #     valid_indices = self.mask[epoch][~np.isnan(self.mask[epoch])].astype(int)
            #     self.mask[epoch] = self.mask[epoch, valid_indices]

            #     mask = self.mask[epoch]
            #     self.pseudo_ranges[epoch][mask] = np.nan
            #     self.elevation[epoch][mask] = np.nan
            # else:
            #     self.mask[epoch] = ~np.isnan(self.mask[epoch])
        for i, obs in enumerate(self.pseudo_ranges):
            # self.epoch = i
            # with open("exemplo.txt", "w") as file_txt:      
                self.__update(obs, i)  
                end = time.perf_counter()        
                print(f"epoca {i} tempo: {end - start} s")

        data = {
            "x": self.obs_time,
            "y": self.pos_calc
        }

        with open("sem_corr_data.pkl", "wb") as f:
            pickle.dump(data, f) 
        
            # print(self.Xk)
        # plt.show()
        
