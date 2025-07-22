import numpy as np
class Elevation():
        
        
    def create_rotation_matrix(self, lat, lon):
        slat = np.sin(lat)
        clat = np.cos(lat)
        slon = np.sin(lon)
        clon = np.cos(lon)
        
        R = np.array([
            [-slon, clon, 0],
            [-slat * clon, -slat * slon, clat],
            [clat * clon, clat * slon, slat]
        ])
        return R

    def __compute_coord_enu(self, receiver_coord, sat_coords):
        dif_sat = []

        
        # Calcula a diferença das coordenadas
        dif_sat = sat_coords - receiver_coord
        
        # Converte as coordenadas do receptor para geodésicas
        lat, lon, _ = self.compute_geodesic_coord(*receiver_coord)
        latitude = np.deg2rad(lat) 
        longitude = np.deg2rad(lon)

        # Cria a matriz de rotação
        R = self.create_rotation_matrix(latitude, longitude)

        # Calcula as coordenadas ENU
        ENU = np.dot(R, dif_sat.T).T

        return ENU, R


    def compute_geodesic_coord(self, x, y, z):
        # Constantes

        major_axis = 6378137.000
        flatt = 1 / 298.257223563
        second_ecc = (2 * flatt) - (flatt ** 2)

        # x, y, z = cartesian_coord

        # Cálculo da longitude
        lon = np.arctan2(y, x)
        lon = np.rad2deg(lon) 

        # Cálculos intermediários
        coord1 = x ** 2
        coord2 = y ** 2
        coord3 = z ** 2
        r = np.sqrt(coord1 + coord2 + coord3)
        p = np.sqrt(coord1 + coord2)
        latg = np.arctan2(z, p)
        lat = latg
        latnow = lat + 0.01
        h = p / np.cos(lat) - major_axis

        # Iteração para calcular a latitude precisa
        while abs(lat - latnow) > 1e-10:
            latnow = lat
            sin2lat = np.sin(latnow) ** 2
            N = major_axis / np.sqrt(1 - (second_ecc * sin2lat))
            h = (p / np.cos(latnow)) - N
            lat = np.arctan((z / p) * (1 / (1 - second_ecc * (N / (N + h)))))

        # Converte latitude para graus
        lat = np.rad2deg(lat)

        # Saída das coordenadas geodésicas
        geodesic_coord = [lat, lon, h]
        return geodesic_coord
    
    def compute_elevation(self, receiver_coord, sat_coords):
        # Computa as coordenadas ENU
        ENU, _ = self.__compute_coord_enu(receiver_coord, sat_coords)
        
        # Calcula as normas
        coord1 = ENU[:, 0] ** 2
        coord2 = ENU[:, 1] ** 2
        coord3 = ENU[:, 2] ** 2
        norm = np.sqrt(coord1 + coord2 + coord3)
        
        # Calcula a Elevação e o Azimute
        Elevation = np.degrees(np.arcsin(ENU[:, 2] / norm))
        Azimuth = np.degrees(np.arctan2(ENU[:, 0], ENU[:, 1]))
        
        # Ajusta o Azimute para o intervalo correto
        Azimuth[Azimuth < 0] += 360
        
        return Azimuth, Elevation
