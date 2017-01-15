# -*- coding: utf-8 -*-

import math
class PVTG:
    def Ppc(self,spGr):
        return 756.8-131.07*spGr-3.6*spGr**2
    def Tpc(self,spGr):
        return 169.2+349.5*spGr-74.0*spGr**2
    def Pr(self,P,Ppc):
        return P/Ppc
    def Tr(self,T,Tpc):
        return T/Tpc


    def gasDensity(self,P,T,spGr):
        ''' pressure in psia, T in Rankine
        returns gas Density in g/cc'''
        Mg = 28.967*spGr
        pseudoReducedPress = self.Pr(P,self.Ppc(spGr))
        pseudoReducedTemp = self.Tr(T,self.Tpc(spGr))
        z = self.zFactor(pseudoReducedPress,pseudoReducedTemp)
        return 0.00149406*P*Mg/(z*T)

    ''' The Lee et al. 20 correlation should be used for gases having a specific gravity less than
    0.77. The ranges of applicability are as follows: 100 ... p (psia) ... 8000 ; 200 ... T(°F) ... 340 ;
    0.55 ... N 2 (mol%) ... 4.8 ; and 0.90 ... CO 2 (mol%) ... 3.20 . The reported accuracy of the Lee
    et al correlation was a standard deviation of  2.69% and a maximum deviation of 8.99%.'''
    '''The Sutton 7 correlation should be used for gases having gravity as high as 1.861'''
    def gasViscosity(self,P,T,spGr,correlation):
        ''' Lee et al. estimation of natural gas viscosity
        P in psia, T in Rankine
        '''
        Mg = 28.967*spGr
        gas_density = self.gasDensity(P,T,spGr)
        if correlation == 'lee':
            K1 = ((10**-4)*(9.379+0.01607*Mg)*T**1.5)/(209+19.26*Mg+T)
            X = 3.448+(986.4/T)+(0.01009*Mg)
            Y = 2.447-0.2224*X
        elif correlation == 'suton':
            X = 3.47+(1588/T)+(0.0009*Mg)
            Y = 1.66378-0.04679*X
            Tc = self.Tpc(spGr)
            Tpr = self.Tr(T,Tc)
            vis_norm=0.9490*((Tc)/(Mg**3)*(self.Pc)**4)**(1/6)
            K1= (((10)**-4)*((.807*(Tpr)**0.618)-(0.357*math.exp(-0.449*Tpr))+(0.340*math.exp(-4.058*Tpr))+0.018))/vis_norm

        return K1*math.exp(X*gas_density**Y)


    def gasFVF(self,P,T,z,flag):
        '''

        flag:
        1 = Bg in rcf/scf
        2 = Bg in RB/scf
        3 = Bg in Rm3/Sm3

        Important:

        if flag=1 or flag=2 then P in psia and T in R
        if flag=3 then P is in Kpa and T in K
        '''
        if flag==1:
            return 0.02819*z*T/P
        elif flag==2:
            return 0.005021*z*T/P
        elif flag==3:
            return 0.350958*z*T/P

    '''The DAK correlation is based on the 11-parameter Starling equation of state. Using non-
    linear regression methods'''
    def zFactor(self,Pr,Tr):
        '''Dranchuk and Abou-Kassem fit of Standing and Katz gas compressibility factor '''
        A1 = 0.3265
        A2 = -1.0700
        A3 = -0.5339
        A4 = 0.01569
        A5 = -0.05165
        A6 = 0.5475
        A7 = -0.7361
        A8 = 0.1844
        A9 = 0.1056
        A10 = 0.6134
        A11 = 0.7210

        zguess=1.0
        delta=1.0
        while (delta>0.0001):
            rhor = 0.27*Pr/(zguess*Tr)
            z = 1+(A1+(A2/Tr)+(A3/Tr**3)+(A4/Tr**4)+(A5/Tr**5))*rhor+((A6+(A7/Tr)+(A8/Tr**2))*rhor**2)-(A9*((A7/Tr)+(A8/Tr**2))*rhor**5)+A10*(1+(A11*rhor**2))*((rhor**2)/(Tr**3))*(math.exp(-A11*(rhor*2)))
            delta=abs(z-zguess)
            zguess=z
        return z


    def rockCompressibility(self,porosity):
        ''' Rock compressibility according to Hall's correlation (1/psi)
         (1953) published a relation between pore volume compressibility,Cpc,and initial porosity,φ'''
        return 1.78e-5*(porosity)**-0.4358

