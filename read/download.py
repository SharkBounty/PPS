from ftplib import FTP
import os
import datetime 
import hatanaka
import sys
import zipfile


class GNSSDataDownloader:
    def __init__(self, station_names):
        self.dirs = []
        self.station_names = station_names
        self.ftp_server = 'gssc.esa.int'
        self.ftp_user = 'anonymous'
        
        self.ftp = FTP(self.ftp_server)
        self.ftp.login(user=self.ftp_user)


    def download_gnss_data(self, d1, d2):
        self.path = f'cddis/gnss/data/daily/{d1[0]}'
        self.check_date_inconsistency(d1[0], d1[1], d1[2])
        self.check_date_inconsistency(d2[0], d2[1], d2[2])

        day_start = self.date_to_doy(d1[0], d1[1], d1[2])
        day_end = self.date_to_doy(d2[0], d2[1], d2[2])

        #check input inconsistencies:

        # self.dirs.append(outPath)

        
            
        self.ftp.cwd(self.path)
        for i in range(day_start, day_end + 1):
            outPath = f"data/{d1[0]}/{i}"
            if not(os.path.exists(outPath)):
                os.makedirs(outPath)
            day = f'{i:03}'

            self.ftp.cwd(day)
            files = self.ftp.nlst()
            matching_files = [file for file in files if (file.startswith('brdc') and (file.endswith('n.Z') or file.endswith('n.gz')))]
            for file in matching_files:
                local_filepath = os.path.join(outPath, file)
                with open(local_filepath, 'wb') as local_file:
                    print(f"Downloading {file} to {local_filepath}...")
                    self.ftp.retrbinary(f"RETR {file}", local_file.write)
            self.ftp.cwd('..')
        self.ftp.quit()

    def download_gnss_data_specific_days(self, days):
        self.ftp.cwd(self.path)
        
        for i in days:
            day = f'{i:03}'
            day_dir = os.path.join('data', day)
            os.makedirs(day_dir, exist_ok=True)
            self.ftp.cwd(day)
            files = self.ftp.nlst()
            for station in self.station_names:
                matching_files = [file for file in files if (file.startswith('brdc') and (file.endswith('n.Z') or file.endswith('n.gz')))]
                for file in matching_files:
                    local_filepath = os.path.join(day_dir, file)
                    with open(local_filepath, 'wb') as local_file:
                        print(f"Downloading {file} to {local_filepath}...")
                        self.ftp.retrbinary(f"RETR {file}", local_file.write)
            self.ftp.cwd('..')
        self.ftp.quit()

    @staticmethod
    def date_to_doy(year, month, day):
        date = datetime.datetime(year, month, day)
        doy = date.timetuple().tm_yday
        return doy

    def unzip(self):
        for root, dirs, files in os.walk('data'):
            for filename in files:
                file_path = os.path.join(root, filename)
                if filename.endswith('.Z') or filename.endswith('.gz'):       
                    try:
                        hatanaka.decompress_on_disk(file_path, delete=True)
                        print(f"Extracted and deleted: {file_path}")
                    except Exception as e:
                        print(f"Failed to extract {file_path}: {e}")
                elif filename.endswith('.zip'):
                    try:
                        with zipfile.ZipFile(file_path, 'r') as zip_ref:
                            filtered_files = [file for file in zip_ref.infolist() if file.filename.endswith(('o', 'm', 'n'))]
                            zip_ref.extract(filtered_files[0], root)
                        print(f"Files extracted to {file_path}")
                        os.remove(file_path)
                    except Exception as e:
                        print(f"Failed to extract {file_path}: {e}")
    def RBMC_list(self):
        """destinada a gerar a lista de todas as estações que já existiram da RBMC"""

        #acessando e realizando login
        ftp = FTP("geoftp.ibge.gov.br")
        ftp.login()

        #acessando o diretorio "relatorio", no qual será possível obter os nomes das estações
        ftp.cwd("informacoes_sobre_posicionamento_geodesico")

        ftp.cwd("rbmc")
        ftp.cwd("relatorio")

        stationslist = []

        ftp.dir(stationslist.append)

        #we can now close the connection
        ftp.quit()

        #removing if not contains a station name
        stationslist[:] = [item for item in stationslist if item.find("escritivo") != -1]

        # only acronyms now
        temp = []
        key = "escritivo"
        for item in stationslist:
            pos = item.find(key)+len(key)
            stationName = item[pos+1:pos+4+1]
            temp.append(stationName)

        stationslist = temp

        #check for inconsistencies
        if stationslist == []:
            sys.exit("lista de estações vazia")

        for item in stationslist:
            if len(item) != 4:
                sys.exit("problemas com a lista de nomes das estações, possível mudança de estrutura do servidor do IBGE, contate os desenvolvedores")

        return stationslist

    def check_date_inconsistency(self, year, month, day):
        """testa se uma data entrada é válida"""

        # errDate = ["entrada inválida para o Dia","entrada inválida para o Mês",
        # "entrada inválida para o Ano",
        # "data negativa fornecida como entrada"]

        # # non integer
        # if not isinstance(self.day,int):
        #     print(errDate[0])
        #     return True

        # elif not isinstance(self.month,int):
        #     print(errDate[1])
        #     return True

        # elif not isinstance(self.year,int):
        #     print(errDate[2])
        #     return True

        # # non positive
        # elif self.day < 0 or self.month < 0 or self.year < 0:
        #     print(errDate[3])
        #     return True
        
        # #impossible dates
        # elif self.day > 31:
        #     print(errDate[0])
        #     return True
        #     #TODO: a month-dependant check
        
        # elif self.month > 12:
        #     print(errDate[1])
        #     return True

        # elif self.year < 2010 or self.year > datetime.date.today().year:
        #     print(errDate[3])
        #     return True
        
        # else:
        #     return False
        try:
            # print([self.day,self.month,self.year])
            theDate = datetime.datetime(year,month,day)
            
            #we don't have future data
            if theDate > datetime.datetime.today():
                return True
            return False
        except:
            return True

    def download_station_day(self, year, month, day, rinex3=False,download_all=False):
        """principal função, todas as demais chamam essa, uma ou mais vezes. 
        Pode ser usada também pra baixar todas as estações disponíveis para um determinado dia"""
        year = year
        month = month
        day = day
        doy = self.date_to_doy(year, month, day)
        #check input inconsistencies:
        outPath = f"data/{year}/{doy}"
        self.dirs.append(outPath)

        os.makedirs(outPath, exist_ok=True)

        if len(self.station_names) != 4 and not download_all:
            print("deve-se entrar com a sigla da estação, consulte no site do IBGE")
            return

        if self.check_date_inconsistency(year, month, day): #para a data
            print("data fornecida inconsistente")
            return

        #retrieve stations list:
        stations = self.RBMC_list()

        station_name = self.station_names.upper() #if the user enters as a lowercase or mixed

        if station_name in stations or download_all:

            filesList = []

            ftp = FTP("geoftp.ibge.gov.br")
            ftp.login()

            ftp.cwd("informacoes_sobre_posicionamento_geodesico")

            ftp.cwd("rbmc")

            #date as a structure:
            theDate = datetime.datetime(year,month,day)

            #Day Of Year, also as string
            DOY = theDate.timetuple().tm_yday
            # print(DOY)

            DOYstr = str(DOY)
            if DOY <= 9:
                DOYstr = "00"+DOYstr
            elif DOY <= 99:
                DOYstr = "0"+DOYstr

            if not rinex3:
                ftp.cwd("dados")
                try:
                    ftp.cwd(str(year))
                    #######################
                    ftp.cwd(DOYstr)
                    filesList = ftp.nlst()
                    downList = [item for item in filesList if station_name in item]
                    if downList == []:
                        downList = [item for item in filesList if station_name.lower() in item]

                    if(download_all):
                        downList = filesList
                    if downList == []:
                        print("estação indisponivel para a data ou previamente desativada")
                        return
                    else:
                        # print(downList)
                        for file in downList:
                            local_filename = os.path.join(outPath,file)
                            try:
                                folder = open(local_filename, "wb")
                                ftp.retrbinary("RETR " + file, folder.write, 8*1024)

                                folder.close()
                            except:
                                print("problema ao tentar acessar o caminho especificado")
                                return

                except:
                    print("ano indisponível")
                    return
            else:
                ftp.cwd("dados_RINEX3")
                try:
                    ftp.cwd(str(year))
                    if year < 2018:
                        print("dia do ano indisponível, para RINEX 3, dados somente a partir do dia 242 de 2018")
                        return
                    elif DOY < 242 and year <= 2018:
                        print("dia do ano indisponível, para RINEX 3, dados somente a partir do dia 242 de 2018")
                        return
                    else:
                        ######################
                        ftp.cwd(DOYstr)
                        filesList = ftp.nlst()

                        # print(filesList)

                        downList = [item for item in filesList if station_name in item]
                        if downList == []:
                            downList = [item for item in filesList if station_name.lower() in item]

                        
                        if(download_all):
                            downList = filesList

                        if downList == []:
                            print("estação indisponivel para a data, previamente desativada ou sem suporte para RINEX3")
                            return
                        else:
                            #######################
                            for file in downList:
                                local_filename = os.path.join(outPath,file)
                                try:
                                    folder = open(local_filename, "wb")
                                    ftp.retrbinary("RETR " + file, folder.write, 8*1024)

                                    folder.close()
                                except:
                                    print("problema ao tentar acessar o caminho especificado")
                                    return

                except:
                    print("ano indisponível, para RINEX 3, dados somente a partir do dia 242 de 2018")
                    return

            ftp.quit()

        else:
            print("nome de estação inválido, entre com as siglas, consulte no site do IBGE")
            return

    def downBetw2dates(self, d1, d2, rinex3=False,download_all=False):
        """faz download para todos os dias no intervalo entre duas datas, da pra usar pra baixar de todas as estações disponíveis pra um dia"""
        
        try:
            date1 = datetime.datetime(d1[0], d1[1], d1[2])
        except:
            date1 = datetime.datetime(d1[0], d1[1], d1[2])
            print("data de início inválida")
            return
        
        try:
            date2 = datetime.datetime(d2[0],d2[1],d2[2])
        except:
            print("data de fim inválida")
            return
                   

        dateList = []

        if date2 > date1:
            interval = date2 - date1
            days = interval.days

            for num in range (0,interval.days+1):
                dateList.append(date1 + datetime.timedelta(days = num))
            
        else:
            interval = date1 - date2
            days = interval.days

            for num in range (0,interval.days+1):
                dateList.append(date2 + datetime.timedelta(days = num))

        for date in dateList:
            self.download_station_day(date.year, date.month, date.day, rinex3, download_all)    

    def download_station_list_date(self, station_list: list,rinex3=False):
        """faz download de uma lista de estações para uma data"""
        for item in station_list:
            self.download_station_day(item, self.outPath,rinex3)
