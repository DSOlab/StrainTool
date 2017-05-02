def parse_ascii_input(filename):
    stations = []
    with open(filename) as fin:
        for line in fin.readlines():
            stations.append(Station(line))
    #print '[DEBUG] Read {} stations.'.format(len(stations))
    if len(stations):
        return stations
    else:
        return None
