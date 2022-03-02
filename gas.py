#gas_flow_class
class Gas():
    def __init__(self) -> None:
        self.Pt = 0
        self.Tt = 0
        self.P = 0
        self.T = 0
        self.H = 0
        self.S = 0
        self.Cp = 1.005
        self.f = 0 #fraction of air and gas
        self.mass = 0 #massflow
        self.c = 0 #velocity of gas in m/s
        self.a = 0 #sound speed
        self.k = 1.33
        self.Rg = 287
        self.psi = 0