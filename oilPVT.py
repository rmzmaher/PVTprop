import sympy
import math
from sympy import ln

class PVTO:
    # step1: Correcting Gas Gravity to separator conditions
    def Gamma_gsp_correction_VB(self,Ps, Gamma_API, Ts, Gamma_g):
        """
              Gamma_gsp_correction_VB(Separator Pressure(Psia), Oil API Gravity, Separator Temperature (F), Separator Gas Gravity) -> Corrected Gas Gravity
              Is used to correct gas gravity to separator conditions
              Vazquez and Beggs Correction
        """
        Gamma_gsp = Gamma_g * (1 + 5.912 * (10 ** -5) * Gamma_API * Ts * ((ln(Ps / 114.7)) / 2.3026))

    # step2: Calculating Bubble point pressure

    def Pb_Valco_McCain(self,Rsb, Gamma_API, T, Gamma_gsp):
        """
             Pb_Valco_McCain(solution GOR at Pb(scf/STB), Oil API Gravity, Temperature (F), Separator Gas Gravity) -> Bubble point Pressure
             Is used to calculate the oil's Bubble Point pressure from field recordings
             Valco McCain Correlation
        """
        z1 = -5.48 - 0.0378 * ln(Rsb) + 0.281 * ((ln(Rsb)) ** 2) - 0.0206 * ((ln(Rsb)) ** 3)
        z2 = 1.27 - 0.0449 * Gamma_API + 4.36 * (10 ** -4) * (Gamma_API ** 2) - 4.76 * (10 ** -6) * (Gamma_API ** 3)
        z3 = 4.51 - 10.84 * Gamma_gsp + 8.39 * (Gamma_gsp ** 2) - 2.34 * (Gamma_gsp ** 3)
        z4 = -0.7835 + 6.23 * (10 ** -3) * T - 1.22 * (10 ** -5) * (T ** 2) + 1.03 * (10 ** -8) * (T ** 3)
        z = z1 + z2 + z3 + z4
        Pb = math.exp(7.475 + 0.713 * z + 0.0075 * (z ** 2))


    def Pb_Standing(self,Rsb, Gamma_API, T, Gamma_gsp):
        """
            Pb_Standing(solution GOR at Pb(scf/STB), Oil API Gravity, Temperature (F), Separator Gas Gravity) -> Bubble point Pressure
            Is used to calculate the oil's Bubble Point pressure from field recordings
            Standing Correlation
        """
        Cpb = ((Rsb / Gamma_gsp) ** 0.83) * (10 ** (0.00091 * T - 0.0125 * Gamma_API))
        Pb = 18.2 * (Cpb - 1.4)

    # step3: Calculating equilibrium GOR
    # step4: Calculating oil compressibility

    def Co_Spivey(self,Pb,P,Rsb, Gamma_API, T, Gamma_gsp):
        """
             Co_Spivey(Bubble point pressure (psia), Pressure (psia), solution GOR (scf/STB),Oil API Gravity, Temperature (F), Separator Gas Gravity) -> oil compressibility
             Is used to calculate the oil's compressibility from field recordings
             Spivey etal. Correlation
        """
        z1 = 3.011-2.6254*ln(Gamma_API)+0.497*((ln(Gamma_API))**2)
        z2 = -0.0835-0.259*ln(Gamma_gsp)+0.382*((ln(Gamma_gsp))**2)
        z3 = 3.51 - 0.0289*ln(Pb)-0.0584*((ln(Pb))**2)
        z4 = 0.327 - 0.608*ln(P/Pb)+0.0911*((ln(P/Pb))**2)
        z5 = -1.918-0.642*ln(Rsb)+0.154*((ln(Rsb))**2)
        z6 = 2.52-2.73*ln(T)+0.426*((ln(T))**2)
        z=z1+z2+z3+z4+z5+z6
        Cob=math.exp(2.434+0.475*z+0.048*(z**2))

    def Co_Marhoun(self,Pb,P,Rsb, Gamma_o, T, Gamma_gsp,Bob):
        """
             Co_Marhoun(Bubble point pressure (psia), Pressure (psia), solution GOR (scf/STB),Oil API Gravity, Temperature (F), Average Gas Gravity) -> oil compressibility
             Is used to calculate the oil's compressibility from field recordings
             El Marhoun Correlation
        """
        Gamma_ob = (Gamma_o+2.18*(10**-4)*Rsb*Gamma_gsp)/Bob
        Co = math.exp(-14.1042+(2.7314/Gamma_ob)-56.0605*(10**-6)*((P-Pb)/(Gamma_ob**3))-580.8778/(T+460))

    # step5: Calculating oil viscosity

    def oil_visc_BR(self,P,Pb,Gamma_API, T, Rs):
        """
               oil_visc_BR(Oil API Gravity, Temperature (F), solution GOR (scf/STB)) -> Oil viscosity (cp)
               is used to calculate the oil's viscosity
               below bubble point:
               The Ng-Egbogah / The Beggs and Robinson correlation
               above bubble point:
               Bergman and Sutton
        """
        Muod = (10 ** ((10 ** (1.8653 - 0.025086 * Gamma_API - 0.5644 * (math.log10(T)))))) - 1
        a = 10.715 * ((Rs + 100) ** -0.515)
        b = 5.44 * ((Rs + 150) ** -0.338)
        Muob = a * (Muod ** b)
        if P <= Pb:
            Muo = Muob
        else:
            a = 6.5698 * (10 ** -7) * ln((Muob) ** 2) - 1.48211 ** (10 ** -5) * ln(Muob) + 2.27887 * (10 ** -4)
            b = 2.24623 * (10 ** -2) * ln(Muob) + 0.873204
            Muo = Muob * math.exp(a * ((P - Pb) ** b))


    # step6: Calculating oil FVF
    def FVFo(self,P,Pb,Gamma_g,Gamma_o, T, Rs):
        if P<=Pb:
            Bo = 1.0113+7.2046*(10**-5)*(((Rs**0.3738)*((Gamma_g**0.2914)/(Gamma_o**0.6265))+0.2462*(T**0.5371))**3.0936)




