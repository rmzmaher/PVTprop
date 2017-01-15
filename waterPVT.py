
''' The water properties required for most well test application
are FVF Bw , Viscosity mu_w and cw'''
class PVTW:

    def waterCompressibility(self,P,T,salinity,correlation):
        ''' Water compressibility in psi**-1 (Osif's laboratory measurements)
        P = pressure, psi
        T = temperature, F
        salinity in mg/l of solution (good for 0 to 200,000 mg/l NaCl)
        '''
        if correlation == 'osif':
            m1 = 7.033
            m2 = 0.5415
            m3 = -537
            m4 = 403300
            cwInv = m1*P+m2*salinity+m3*T+m4
            return 1/cwInv

    def waterViscosity(self,P,T,spGr,correlation):
        ''' McCain: estimation of water viscosity
        P in psia, T in Rankine
        '''

    def waterFVF(self,P,T,correlation):
        '''
        McCain correlation :The correlation is valid for temperatures up to 260 ° F and pressures up to 5000 psia
        :return:
        '''
        if correlation == 'McCain':
            DVwT = - 1.0001 * (10**-2) + 1.33391 * (10**-4) * T + 5.50654 * (10 ** -7) * T**2
            DVwp = - 1.95301 * (10**-9) * P * T - 1.72834 * 10**-13 * (P**2) * T - 3.58922 * (10**-7) * P - 2.25341 * (10**-10) * (P**2)
            Bw = (1 + DVwp )(1 + DVwT )
            return Bw
