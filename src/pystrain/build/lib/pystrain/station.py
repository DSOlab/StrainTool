class Station:
    def __init__(self, input_line):
        ''' parse from standard ascii input line '''
        l = input_line.split()
        self.name = l[0]
        self.lon  = float(l[1])
        self.lat  = float(l[2])
        self.ve   = float(l[3])
        self.vn   = float(l[4])
        self.se   = float(l[5])
        self.sn   = float(l[6])
        self.rho  = float(l[7])
        self.t    = float(l[8])
