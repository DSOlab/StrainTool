from pystrain.station import Station

def parse_ascii_input(filename, zero_std_is_error=False):
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
          zero_std_is_error (bool): if set to True, then the function will throw
                             if a station has zero std. deviation for either the
                             north or east component (or both)

      Returns:
          list of Station instances or None (if no station was read)

  """
  stations = []
  with open(filename) as fin:
    for line in fin.readlines():
      # stations.append(Station(line))
      nSta=Station(line)
      if zero_std_is_error and (nSta.sn==0e0 or nSta.se==0e0):
        raise ValueError('[ERROR] Zero std. deviation not allowed! station is: {:}'.format(nSta.name))
      ## check that the station is not a duplicate
      for sta in stations:
        if sta.name == nSta.name:
          raise ValueError('[ERROR] Duplicate record found in input file for station {:}'.format(sta.name))
        if sta.lat==nSta.lat and sta.lon==nSta.lon:
          raise ValueError('[ERROR] Exact coordinate match for stations {:} and {:}. Possible duplicate!'.format(sta.name, nSta.name))
      stations.append(nSta)
  if len(stations):
    return stations
  else:
    return None
