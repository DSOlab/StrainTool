from pystrain.station import Station

def parse_ascii_input(filename):
    stations = []
    with open(filename) as fin:
        for line in fin.readlines():
            stations.append(Station(line))
    if len(stations):
        return stations
    else:
        return None
