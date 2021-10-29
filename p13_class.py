import numpy as np

def FNday(Y, M, D):
    YM  = 365.25             # Mean Year,     days
    if M <= 2:
        Y = Y - 1
        M = M + 12
    return int(Y * YM) + int((M + 1) * 30.6) + D - 428

class plan13:
    def __init__(self):
        SAT = "OSCAR-13"
        YE      = 1990           # Epoch Year
        self.TE = 191.145409     # Epoch time (days)
        IN      =  56.9975       # Inclination   deg
        RA      = 146.4527       # R.A.A.N.      deg
        self.EC =   0.6986       # Eccentricity   -
        WP      = 231.0027       # Arg perigee   deg
        MA      =  43.2637       # Mean anomaly  deg
        MM      =   2.09695848   # Mean motion   rev/d
        M2      = 1E-8           # Decay Rate    rev/d/d
        self.RV = 1585           # Orbit number   -
        ALON    = 180            # Sat attitude, deg. 180 = nominal ) See bulletins
        ALAT    =   0            # Sat attitude, deg.   0 = nominal ) for latest

        # Observer's location + North, + East, ASL(m) (Geodetic)
        LOC="G3RUH"
        LA = np.deg2rad(52.21)
        LO = np.deg2rad(0.06)
        HT = 79 / 1000.0

        """ *******************************************************************************************
        Converting geodetic to ECEF, see https://en.wikipedia.org/wiki/Geographic_coordinate_conversion
        ******************************************************************************************* """
        
        CL = np.cos(LA)
        SL = np.sin(LA)
        CO = np.cos(LO)
        SO = np.sin(LO)

        # WGS-84 Earth ellipsoid
        self.RE = 6378.137
        FL = 1/298.257224

        RP = self.RE*(1-FL)
        XX = self.RE*self.RE
        ZZ = RP*RP

        D  = np.sqrt(XX*CL*CL + ZZ*SL*SL)
        Rx = XX/D + HT
        Rz = ZZ/D + HT

        # Observer's unit vectors UP EAST and NORTH in GEOCENTRIC coords.
        # Up XYZ
        self.Ux = CL*CO
        self.Uy = CL*SO
        self.Uz = SL

        # East XYZ
        self.Ex = -SO
        self.Ey = CO
        self.Ez =  0

        # North XYZ
        self.Nx = -SL*CO
        self.Ny = -SL*SO
        self.Nz = CL

        # Observer's XYZ coords at Earth's surface
        self.Ox = Rx*self.Ux
        self.Oy = Rx*self.Uy
        self.Oz = Rz*self.Uz

        # Convert angles to radians etc.
        self.RA = np.deg2rad(RA)
        IN = np.deg2rad(IN)
        self.WP = np.deg2rad(WP)
        self.MA = np.deg2rad(MA)
        self.MM = MM*2*np.pi
        M2 = M2*2*np.pi

        YM = 365.25             # Mean Year,     days
        YT = 365.2421874        # Tropical year, days
        self.WW = 2*np.pi/YT       # Earth's rotation rate, rads/whole day
        self.WE = 2*np.pi + self.WW     #       ditto            radians/day
        W0 = self.WE/86400           #       ditto            radians/sec

        # Observer's velocity, GEOCENTRIC coords. (VOz=0)
        self.VOx = -self.Oy*W0
        self.VOy = self.Ox*W0

        # Convert satellite Epoch to Day No. and Fraction of day
        self.DE = self.FNday(YE, 1, 0) + int(self.TE)
        self.TE = self.TE - int(self.TE)

        """ Average Precession rates """
        GM      = 3.986E5                # Earth's Gravitational constant km^3/s^2
        J2      = 1.08263E-3             # 2nd Zonal coeff, Earth's Gravity Field
        self.N0 = self.MM/86400               # Mean motion rad/s
        self.A0 = (GM/self.N0/self.N0)**(1/3)       # Semi major axis km
        self.B0 = self.A0*np.sqrt(1-self.EC*self.EC)  # Semi minor axis km
        self.SI = np.sin(IN)
        self.CI = np.cos(IN)
        PC      = self.RE*self.A0/(self.B0*self.B0)
        PC      = 1.5*J2*PC*PC*self.MM        # Precession const, rad/Day
        self.QD = -PC*self.CI                 # Node precession rate, rad/day
        self.WD =  PC*(5*self.CI*self.CI-1)/2      # Perigee precession rate, rad/day
        self.DC = -2*M2/self.MM/3             # Drag coeff. (Angular momentum rate)/(Ang mom)  s^-1

        """ Please see end of listing for newer values; use old ones for test. """
        # Sidereal and Solar data. Rarely needs changing. Valid to year ~2015
        YG = 2000
        G0 = 98.9821 # GHAA, Year YG, Jan 0.0
        # MA Sun and rate, deg, deg/day
        MAS0 = 356.0507
        self.MASD = 0.98560028
        # Sun's inclination
        INS = np.radians(23.4393)
        self.CNS = np.cos(INS)
        self.SNS = np.sin(INS)
        # Sun's Equation of centre terms
        self.EQC1 = 0.03342
        self.EQC2 = 0.00035

        # Bring Sun data to Satellite Epoch
        TEG  = (self.DE - self.FNday(YG,1,0)) + self.TE            # Elapsed Time: Epoch - YG
        self.GHAE = np.radians(G0) + TEG*self.WE            # GHA Aries, epoch
        self.MRSE = np.radians(G0) + TEG*self.WW + np.pi  # Mean RA Sun at Sat epoch
        self.MASE = np.radians(MAS0 + self.MASD*TEG)        # Mean MA Sun  ..

        # Antenna unit vector in orbit plane coordinates.
        CO = np.cos(np.radians(ALON))
        SO = np.sin(np.radians(ALON))
        CL = np.cos(np.radians(ALAT))
        SL = np.sin(np.radians(ALAT))
        self.ax = -CL*CO
        self.ay = -CL*SO
        self.az = -SL

        # Miscellaneous
        OLDRN = -99999

    def PROCsatvec(self, T):
        # Linear drag terms
        DT  = self.DC * T/2.0
        KD  = 1.0 + 4.0*DT
        KDP = 1.0 - 7.0*DT

        # Mean anomaly at YR,TN
        M  = self.MA + self.MM*T*(1.0 - 3.0*DT)

        # Strip out whole no of revs
        DR = int(M/(2.0*np.pi))

        # M now in range 0 - 2pi
        M  = M - DR*2.0*np.pi

        # Current Orbit number
        self.RN = self.RV + DR

        """ Solve M = EA - EC*SIN(EA) for EA given M, by Newton's Method """
        EA = M

        while True:
            C = np.cos(EA)
            S = np.sin(EA)
            DNOM = 1 - self.EC*C
            D = (EA - self.EC*S - M)/DNOM   # Change to EA for better solution
            EA = EA - D                # by this amount

            if np.abs(D) < 1E-5:
                break

        # Distances
        A = self.A0*KD
        B = self.B0*KD
        self.RS = A*DNOM

        # Calc satellite position & velocity in plane of ellipse
        Sx = A*(C - self.EC)
        Vx = -A*S/DNOM*self.N0
        Sy = B*S
        Vy = B*C/DNOM*self.N0

        AP = self.WP + self.WD*self.T*KDP
        CW = np.cos(AP)
        SW = np.sin(AP)
        RAAN = self.RA + self.QD*self.T*KDP
        CQ = np.cos(RAAN)
        SQ = np.sin(RAAN)

        # Plane -> celestial coordinate transformation, [C] = [RAAN]*[IN]*[AP]
        CXx =  CW*CQ - SW*self.CI*SQ
        CXy = -SW*CQ - CW*self.CI*SQ
        CXz =  self.SI*SQ

        CYx =  CW*SQ + SW*self.CI*CQ
        CYy = -SW*SQ + CW*self.CI*CQ
        CYz = -self.SI*CQ

        CZx =  SW*self.SI
        CZy =  CW*self.SI
        CZz =  self.CI

        """ Compute SATellite's position vector, ANTenna axis unit vector
            and VELocity in CELESTIAL coordinates. (Note: Sz=0, Vz=0) """
        self.SATx = Sx*CXx + Sy*CXy
        self.SATy = Sx*CYx + Sy*CYy
        self.SATz = Sx*CZx + Sy*CZy

        self.ANTx = self.ax*CXx + self.ay*CXy + self.az*CXz
        self.ANTy = self.ax*CYx + self.ay*CYy + self.az*CYz
        self.ANTz = self.ax*CZx + self.ay*CZy + self.az*CZz

        VELx = Vx*CXx + Vy*CXy
        VELy = Vx*CYx + Vy*CYy
        VELz = Vx*CZx + Vy*CZy

        # Also express SAT,ANT and VEL in GEOCENTRIC coordinates:
        GHAA = self.GHAE + self.WE*T           # GHA Aries at elapsed time T
        C = np.cos(-GHAA)
        S = np.sin(-GHAA)
        
        self.Sx = self.SATx*C - self.SATy*S
        self.Sy = self.SATx*S + self.SATy*C
        self.Sz = self.SATz

        self.Ax = self.ANTx*C - self.ANTy*S
        self.Ay = self.ANTx*S + self.ANTy*C
        self.Az = self.ANTz

        self.Vx = VELx*C - VELy*S
        self.Vy = VELx*S + VELy*C
        self.Vz = VELz

    def PROCsunvec(self):
        MAS = self.MASE + np.radians(self.MASD * self.T)   # MA of Sun round its orbit
        TAS = self.MRSE + self.WW*self.T + self.EQC1*np.sin(MAS) + self.EQC2*np.sin(2*MAS)

        # Sin/Cos Sun's true anomaly
        C = np.cos(TAS)
        S = np.sin(TAS)

        # Sun unit vector - CELESTIAL coords
        self.SUNx = C
        self.SUNy = S*self.CNS
        self.SUNz = S*self.SNS

    def PROCrangevec(self):
        """ Compute and manipulate range/velocity/antenna vectors """
        # Rangevec = Satvec - Obsvec
        Rx = self.Sx - self.Ox
        Ry = self.Sy - self.Oy
        Rz = self.Sz - self.Oz 

        R = np.sqrt(Rx*Rx+Ry*Ry+Rz*Rz)         # Range magnitude
        self.R = R
        # Normalise Range vector
        Rx=Rx/R
        Ry=Ry/R
        Rz=Rz/R 

        U = Rx*self.Ux + Ry*self.Uy + Rz*self.Uz    # UP    Component of unit range
        E = Rx*self.Ex + Ry*self.Ey                 # EAST    do   (Ez=0)
        N = Rx*self.Nx + Ry*self.Ny + Rz*self.Nz    # NORTH   do
        AZ = np.rad2deg(np.arctan2(E, N))           # Azimuth
        EL = np.rad2deg(np.arcsin(U))               # Elevation
        self.EL = EL
        self.AZ = AZ

        # Resolve antenna vector along unit range vector, -r.a = Cos(SQ)
        SQ = np.rad2deg(np.arccos(-(self.Ax*Rx + self.Ay*Ry + self.Az*Rz))) # Hi-gain ant SQuint
        self.SQ = SQ

        # Calculate sub-satellite Lat/Lon
        SLON = np.rad2deg(np.arctan2(self.Sy, self.Sx))        # Lon, + East
        SLAT = np.rad2deg(np.arcsin(self.Sz/self.RS))    # Lat, + North
        self.SLON = SLON
        self.SLAT = SLAT

        # Resolve Sat-Obs velocity vector along unit range vector. (VOz=0)
        RR  = (self.Vx - self.VOx)*Rx + (self.Vy - self.VOy)*Ry + self.Vz*Rz  # Range rate, km/s
        self.RR = RR

    def isSatIlluminated(self):
        """ Find Solar angle, illumination, and eclipse status. """
        SSA = -(self.ANTx*self.SUNx + self.ANTy*self.SUNy + self.ANTz*self.SUNz)            # Sin of Sun angle -a.h
        ILL = np.sqrt(1 - SSA*SSA)                                                          # Illumination
        CUA = -(self.SATx*self.SUNx + self.SATy*self.SUNy + self.SATz*self.SUNz)/self.RS    # Cos of umbral angle -h.s
        UMD = self.RS*np.sqrt(1 - CUA*CUA)/self.RE                                          # Umbral dist, Earth radii

        if CUA >= 0:
            ECL = "    +"
        else:
            ECL = "    -"

        if (UMD <= 1) and (CUA >= 0):
            ECL = "   ECL"

    def FNday(self, Y, M, D):
        YM  = 365.25             # Mean Year,     days
        if M <= 2:
            Y = Y - 1
            M = M + 12
        return int(Y * YM) + int((M + 1) * 30.6) + D - 428

    def calcSat(self, DN, TN):
        T = (DN - self.DE) + (TN - self.TE) # Elapsed time since epoch, days
        self.T = T
        self.PROCsatvec(T)
        self.PROCsunvec()
        self.isSatIlluminated()
        self.PROCrangevec()


if __name__ == "__main__":
    test = plan13()

    test.calcSat(FNday(1990, 11, 3), 1/24)

    print(test.RN)
    print('0100  {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}'. format(test.R, test.EL, test.AZ, test.SQ, test.RR, test.RS - test.RE, test.SLAT, test.SLON))

    test.calcSat(FNday(1990, 11, 3), 1.25/24)

    print(test.RN)
    print('0115  {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}'. format(test.R, test.EL, test.AZ, test.SQ, test.RR, test.RS - test.RE, test.SLAT, test.SLON))

    test.calcSat(FNday(1990, 11, 3), 1.50/24)

    print(test.RN)
    print('0130  {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}'. format(test.R, test.EL, test.AZ, test.SQ, test.RR, test.RS - test.RE, test.SLAT, test.SLON))

    test.calcSat(FNday(1990, 11, 3), 1.75/24)

    print(test.RN)
    print('0145  {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}'. format(test.R, test.EL, test.AZ, test.SQ, test.RR, test.RS - test.RE, test.SLAT, test.SLON))

    test.calcSat(FNday(1990, 11, 3), 2/24)

    print(test.RN)
    print('0200  {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}'. format(test.R, test.EL, test.AZ, test.SQ, test.RR, test.RS - test.RE, test.SLAT, test.SLON))