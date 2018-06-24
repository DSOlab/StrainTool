from pystrain.station import Station

def parse_ascii_input(filename):
    """Parse station info from an input file.

        This function will try to read Stations of from the input file, named
        filename. For this to work, the lines in the input file, should conform
        to the following format:
        "name lon lat Ve Vn Se Sn RHO T"
        where lon and lat are in decimal degrees and velocity components and
        sigmas are in mm/year.
        For more information, see the functions:
        Station.init_from_ascii_line()
        Station.Station() --aka the constructor--.

        Args:
            filename (string): the name of the file holding station info (see
                               description above).

        Returns:
            list of Station instances or None (if no station was read)

    """
    stations = []
    with open(filename) as fin:
        for line in fin.readlines():
            stations.append(Station(line))
    if len(stations):
        return stations
    else:
        return None
