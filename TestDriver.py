# -*- coding: utf-8 -*-

from RTCoef import RTCoef
import numpy as np
def main():
    """
    usage: In command line, enter "python3 TestDriver.py". And these two files must be in the same directory, if you didn't set the environment variables.

    These medir paremeters like rho(density) and beta(shear wave velocity) is from PREM standrad modle in the depth of 15km.
    """
    icang = [np.deg2rad(i) for i in range(91)]
    """
    beta2=3.2 rho2=2.6
    -----------------------------------------------------------
    beta1=3.9 rho1=2.9  
    
    In addition, R and T coefficient calculated here is not absolute value!! 
    """
    a = RTCoef(*icang, beta1=3.9, beta2=3.2, rho1=2.9, rho2=2.6)
    a.plot()

    """
    beta1=3.2 rho1=2.6
    -----------------------------------------------------------
    beta2=3.9 rho2=2.9
    """

    b = RTCoef(*icang, beta1=3.2, beta2=3.9, rho1=2.6, rho2=2.9)
    b.plot()


if __name__=='__main__':
    main()

