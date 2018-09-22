# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

class RTCoef():
    def __init__(self, *icang, **media):
        """
        icang                   list of incident angle(unit is rad)
        **media(key word argument)
            rho1, rho2          density
            beta1, beta2        S wave velocity

            subscript "1" means the incident side, and "2" means the transmission side.
            
            rho1, beta1
            -----------------------------------
            rho2, beta2

        This class can only calculate the case of solid-solid interface and SH wave.
        """
        self.icang = icang
        self.rho1 = media['rho1']
        self.rho2 = media['rho2']
        self.beta1 = media['beta1']
        self.beta2 = media['beta2']
        if self.beta2 < self.beta1:
            self.isreal = True 
            self.ic = None # ic means critical incident angle.
        else:
            self.isreal = False
            self.ic = np.arcsin(self.beta1 / self.beta2)

        self.RCoef, self.TCoef = self.__calRTCoef()
        #self.RCoef = np.abs(self.RCoef)
        self.RC_Real = [i.real for i in self.RCoef]
        self.RC_Imag = [i.imag for i in self.RCoef]
        self.TC_Real = [i.real for i in self.TCoef]
        self.TC_Imag = [i.imag for i in self.TCoef]
    
    def __calRTCoef(self):
        """
        calculate RTCoef, and the formulas used here is from "Aki and Richards, Quantitative Seismolgy, 2002".
        """
        RC = []
        TC = []
        if self.isreal:
           for i in self.icang:
               tmp1 = self.rho1 * self.beta1 * np.cos(i)
               tmp2 = self.rho2 * self.beta2 * np.cos(self.__gettang(i))
               tmp3 = 2 * tmp1
               RC.append((tmp1 - tmp2) / (tmp1 + tmp2))
               TC.append(tmp3 / (tmp1 + tmp2))
        else:
            for i in self.icang:
                """
                if i < self.ic:
                    tmp1 = self.rho1 * self.beta1 * np.cos(i)
                    tmp2 = self.rho2 * self.beta2 * np.cos(self.__gettang(i))
                    tmp3 = 2 * tmp1
                    RC.append(complex((tmp1 - tmp2) / (tmp1 + tmp2)))
                    TC.append(complex(tmp3 / (tmp1 + tmp2)))
                else:
                    """
                if self.beta2 * np.sin(i) / self.beta1 <= 1:
                    cosj2 = np.sqrt(1 - np.sin(self.__gettang(i)) ** 2)
                else:
                    cosj2 = complex(0, np.sqrt(self.__gettang(i) ** 2 - 1))
                tmp1 = self.rho1 * self.beta1 * np.cos(i)
                tmp2 = self.rho2 * self.beta2 * cosj2
                tmp3 = 2 * tmp1
                RC.append((tmp1 - tmp2) / (tmp1 + tmp2))
                TC.append(tmp3 / (tmp1 + tmp2))
                
        return RC, TC

    def __gettang(self, iang):
        """
        get transmission angle
        """
        if self.isreal:
            return np.arcsin((self.beta2 / self.beta1) * np.sin(iang))
        else:
            if iang < self.ic:
                return np.arcsin((self.beta2 / self.beta1) * np.sin(iang))
            else:
                return np.pi / 2                  


    def plot(self):
        fig = plt.figure(figsize=(10, 4))
        #plt.subplot(211)
        plt.plot(np.rad2deg(self.icang), self.RC_Real, color='red', label="Reflect_Real")
        plt.plot(np.rad2deg(self.icang), self.RC_Imag, color='red', linestyle="--", label="Reflect_Imag")
        #plt.xlabel("Incident Angle")
        #plt.ylabel("R/T coefficient")
        #plt.legend()
        #plt.subplot(212)
        plt.plot(np.rad2deg(self.icang), self.TC_Real, color='blue', label="Transmission_Real")
        plt.plot(np.rad2deg(self.icang), self.TC_Imag, color='blue', linestyle="--", label="Transmission_Imag")        
        plt.legend()
        plt.xlabel("Incident Angle")
        plt.ylabel("R/T coefficient")
        plt.show()
        

