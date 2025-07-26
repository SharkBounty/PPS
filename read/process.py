import numpy as np
import datetime
import xarray as xr


import numpy as np
class Process():
    # def __init__(self, ephemerides):
    #     """
    #     eph    ->  Array contendo as efemérides de uma observável\n 
    #     """
    #     self.data = ephemerides

    def clear_nan_values(self, data):
        if data.ndim > 1:
            available = ~np.isnan(np.array(data)).any(axis=1)
            avail = data[available]  
        else:
            available = ~np.isnan(np.array(data))
            avail = data[available]
        return avail
    
    def get_data_avail(self, data, prefix='G', filter=False, sv_health=False):
        """
        prefix -> Inicial do código do satélite alvo ('E', 'G', 'R')
        """
        
        data_avail = []
        if filter: data = self.filter(data, sv_health)
        for obs in data:
            # available = ~np.isnan(np.array(obs))
            # avail = obs[available]
            target = np.strings.startswith(obs.sv.data, prefix)
            data_avail.append(obs[target])
                
        # merged = xr.concat(data_avail, dim="time")
        return data_avail
    

    def __get_index(self, obs_time, dtime):
        index = []
        for time in obs_time:
            index.append(np.argmin(np.abs(time.data - dtime)))
        return index
            

    def get_data_obs_for_intersect_nav(self, obs, dtime, nav):
        """
        prefix -> Inicial do código do satélite alvo ('E', 'G', 'R')
        """
        
        # sats = [sat.sv.data for sat in obs]
        # times = np.array([sat.sv.time.values for sat in obs])
        obs = xr.concat(obs, dim='time')

        index = self.__get_index(obs.time, dtime)
        
        
        sv_avail = [ list( set(nav[index].sv.data) & set(obs[i].sv.data)) for i, index in enumerate(index)]
        print(sv_avail)
        return sv_avail
    


    def get_tow(self, year1, month1, day1, year2, month2, day2):
        """
        Retorna o tempo em segundos da semana \n
        tow = (D * 24 + H) * 3600 + M * 60 + S \n
        D refere-se ao número do dia da semana – iniciando-se em 0 para o domingo \n
        H, M e S referem-se respectivamente a hora, minuto e segundo da época considerada
        """
        dtime_start = datetime.datetime(year1, month1, day1)
        dtime_end = datetime.datetime(year2, month2, day2)

        days = (dtime_end - dtime_start).days + 1

        dow = dtime_start.toordinal() % 7
        tow = dow*86400
        series = np.arange(tow, tow + 15 * 5760 * days, 15)

        return series

    def date_to_doy(self, year, month, day):
        date = datetime.datetime(year, month, day)
        doy = date.timetuple().tm_yday
        return doy


    def filter(self, data, sv_health):
        var_filtered = []
        prev_avail = []
        for i, var in enumerate(data):
            if not(var.time.dt.minute.item() != 0 and var.time.dt.minute.item() != 59):
                # CORREÇÃO: Acessa sv_health[i] diretamente (é um array numpy)
                # e usa indexação booleana direta, como sugerido pelo erro.
                # Um satélite com saúde 63 é considerado inutilizável.
                is_unhealthy = (sv_health[i] == 63)
                var[is_unhealthy] = np.nan

                if var.time.dt.minute.item() == 59:
                    prev_avail = ~np.isnan(var)
                elif len(prev_avail) > 0:
                    index_avail = np.where((np.isnan(var) & prev_avail))
                    var.data[index_avail] = data[i-1].data[index_avail]
                    # Aplica a máscara de saúde novamente após preencher os dados
                    var[is_unhealthy] = np.nan
                    var_filtered.append(var)
            if i == 0:
                var_filtered.append(var)    
                
        merged = xr.concat(var_filtered, dim="time")
        return merged

    def get_dtime(self, data):
        times = []
        for i in data:
            times.append(i.time.data)
        return times

    # def get_coef_dts(self, nav, obs):
    #     """
    #     nav -> Array contendo os dados de navegação
    #     """
    #     # min_obs_delta_time = 0 #diferença de tempo entre a obtenção de navegação e observação

    #     a0 = []
    #     a1 = []
    #     a2 = []
    #     dtime = []

    #     bias_filtered = self.filter(nav['SVclockBias'])
    #     drift_filtered = self.filter(nav['SVclockDrift'])
    #     driftrate_filtered = self.filter(nav['SVclockDriftRate'])

    #     for i, data in enumerate(bias_filtered):
    #             a0.append(data)
    #             a1.append(drift_filtered[i])
    #             a2.append(driftrate_filtered[i])
    #             dtime.append(data.time.values)

    #     a0 = xr.concat(a0, dim='time', coords='different')
    #     a1 = xr.concat(a1, dim='time', coords='different')
    #     a2 = xr.concat(a2, dim='time', coords='different')
        
    #     coef = xr.Dataset({"a0":a0, 
    #                 "a1":a1,
    #                 "a2":a2,
    #                 "dtime": (('time',), dtime)
    #                 })
    #     return coef

        
    