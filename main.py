
from process.single_point_position import SPP 
from read.download import GNSSDataDownloader
import cProfile
import pyinstrument


if __name__ == '__main__':
   
   
    station = "chpi"

    year = 2014
    month = 9
    day = 27

    start_day = day
    end_day   = start_day



    # downloader = GNSSDataDownloader(station)
    # downloader.download_gnss_data([2014, 9, 26], [2014, 9, 27])

    # downloader.downBetw2dates([2014, 9, 26], [2014, 9, 27])
    # downloader.unzip()

    
    spp = SPP(paths=["data/2014/269", "data/2014/270"], start=[2014, 9, 27], end=[2014, 9, 27], type_obs="C1") 
    # cProfile.run('spp.run()')

    # with pyinstrument.profile():
    spp.run()



