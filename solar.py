from numpy import pi, sin, cos, tan, arccos, sign
from datetime import datetime, timedelta
from pandas import Timestamp

# considerando petrolina
LAT = (-40.507778)
STD_LAT = (-45)


def solartime(stdtime, Lloc = LAT, Lstd = STD_LAT):
    """
    retorna o horário solar como datetime
    stdtime = hora de relógio
    Lloc = longitude local/geográfica em graus
    Lstd = longitude padrão/de fuso horário em graus
    """
    if type(stdtime) not in (datetime, Timestamp):
        raise TypeError()

    n = (stdtime - datetime(stdtime.year, 1, 1)).days
    B = (n - 1) * (2 * pi / 365)
    E = 229.2 * ( 0.000075 + 0.001868 * cos(B) - 0.032077 * sin(B) - 0.014615 * cos(2 * B) - 0.04089 * sin(B))

    return stdtime + timedelta(minutes = 4 * (Lstd - Lloc) + E)



def angHora(stdtime, lat = LAT):
    """
    retorna o ângulo horário em graus no momento
    lat = latitude em graus
    stdtime = momento como timestamp
    """
    if not -90 <= lat <= 90:
        raise ValueError(lat)
    if type(stdtime) not in (datetime, Timestamp):
        raise TypeError(type(stdtime))

    n = (stdtime - datetime(stdtime.year, 1, 1)).days + 1

    nHora = ws(n, lat) / 15
    meioDia = datetime(stdtime.year, stdtime.month, stdtime.day, 12, 0)
    meioDiaSolar = solartime(meioDia)
    horaPassada = (stdtime - meioDiaSolar).seconds / 3600

    return  ws(n, lat) * (horaPassada / nHora)



def declin(n):
    """
    retorna a declinação solar em graus
    n = dia juliano
    """
    if not 0 < n <= 366:
        raise ValueError(n)

    return 23.45 * sin(2 * pi * (284 + n) / 365)



def ws(n, lat = LAT):
    """
    retorna o ângulo horário poente em graus
    n = dia juliano
    lat = latitude em graus
    """
    if not 0 < n <= 366:
        raise ValueError(n)
    if not -90 <= lat <= 90:
        raise ValueError(lat)

    dec = declin(n) * (pi / 180)
    lat = lat * (pi / 180)

    return arccos( -tan(lat) * tan(dec) ) * (180 / pi)



def azimute(stdtime, lat = LAT):
    """
    retorna o azimute solar em graus
    stdtime = momento como timestamp
    lat = latitude em graus
    """
    if not -90 <= lat <= 90:
        raise ValueError(lat)
    if type(stdtime) not in (datetime, Timestamp):
        raise TypeError(type(stdtime))

    n = (stdtime - datetime(stdtime.year, 1, 1)).days + 1

    thetaz = arccos(cos_theta(stdtime, lat))
    w = angHora(stdtime, lat) * (pi / 180)
    dec = declin(n) * (pi / 180)

    return sign(w) * abs(arccos((cos(thetaz) * sin(lat) - sin(dec))/sin(thetaz) * cos(lat))) * (180 / pi)



def cos_theta(stdtime, lat = LAT, beta = 0):
    """
    retorna o coseno de theta
    stdtime = momento como timestamp
    lat = latitude em graus
    beta = inclinação em graus
    """
    if not -90 <= lat <= 90:
        raise ValueError(lat)
    if type(stdtime) not in (datetime, Timestamp):
        raise TypeError(type(stdtime))

    n = (stdtime - datetime(stdtime.year, 1, 1)).days + 1
    w = angHora(stdtime, lat)

    dec = declin(n) * (pi / 180)
    beta = beta * (pi / 180)
    lat = lat * (pi / 180)
    w = w * (pi / 180)

    if beta == 0:
        return cos(lat) * cos(dec) * cos(w) + sin(lat) * sin(dec)
    else:
        az = azimute(stdtime, lat)
        az = az * (pi / 180)
        return sin(dec) * sin(lat) * cos(beta) \
                  - sin(dec) * cos(lat) * sin(beta) * cos(az) \
                   + cos(dec) * cos(lat) * cos(beta) * cos(w) \
                   + cos(dec) * sin(lat) * sin(beta) * cos(az) * cos(w) \
                   + cos(dec) * sin(beta) * sin(az) * sin(w)
