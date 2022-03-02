"""
2021-11-19
此代码提供热力学参数计算，包括：
给定 油气比f 大气湿度d 温度T 计算 比焓H
给定 油气比f 大气湿度d 温度T 计算 熵S
给定 油气比f 大气湿度d 温度T 计算 气体绝热指数k
给定 油气比f 大气湿度d 温度T 计算 热力学常数R
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate


class GasTable():
    def __init__(self):
        self.a_air_i = [-0.30591266e+6, 
               0.10489652e+7, 
               -0.23284057e+6, 
               0.45288431e+6,
               -0.31308477e+6, 
               0.11341362e+6, 
               -0.21298087e+5, 
               0.16363600e+4, 
               0, 
               0.80558643e+4]
        self.a_theta_i = [0, 
               -0.31020306e+6, 
               +0.29961197e+7, 
               -0.27934788e+7, 
               +0.18746407e+7, 
               -0.73499597e+6, 
               +0.15062602e+6, 
               -0.12510984e+5, 
               0, 
               -0.17800633e+4]
        self.a_steam_i = [-0.25179888e+6, 
                 +0.33930851e+5, 
                 +0.43630281e+4, 
                 +0.10594936e+5, 
                 -0.54952291e+4, 
                 +0.14949330e+4, 
                 +0.22772255e+3, 
                 0.17874978e+2, 
                 -0.53236046e+0, 
                 +0.23113968e+3]
        self.__temperature_list= [200, 230, 260, 290, 320, 350, 380, 410, 440, 470, 500, 530, 560, 590, 620, 650, 680, 710, 740, 770, 800, 830, 860, 890, 920, 950, 980, 1010, 1040, 1070, 1100, 1130, 1160, 1190, 1220, 1250, 1280, 1310, 1340, 1370, 1400, 1460, 1490, 1520, 1550, 1580, 1610, 1640, 1670, 1700, 1730, 1760, 1790, 1820, 1850, 1880, 1910, 1940, 2000]
        self.__H_list = [199.92, 229.97, 260.01, 290.06, 320.19, 350.37, 380.63, 410.97, 441.48, 472.07, 502.87, 533.8, 564.93, 596.28, 627.78, 659.56, 691.53, 723.71, 756.1, 788.75, 821.6, 854.62, 887.89, 921.37, 955.02, 988.82, 1022.86, 1057.05, 1091.41, 1125.93, 1160.63, 1195.45, 1230.39, 1265.5, 1300.74, 1336.14, 1371.63, 1407.25, 1442.95, 1478.81, 1514.76, 1586.91, 1623.15, 1659.48, 1695.93, 1732.42, 1769, 1805.66, 1842.4, 1879.23, 1916.1, 1953.06, 1990.09, 2027.17, 2064.33, 2101.54, 2138.83, 2176.2, 2251.07]
        self.__theta_h_list = [69.1, 85.8, 106.3, 129.7, 156.1, 185, 215.9, 240, 284.2, 321, 359.1, 398.8, 439.8, 482.1, 525.6, 570, 615.6, 662.1, 709.8, 758.3, 808.3, 859.6, 911.5, 964.2, 1018.2, 1073, 1129.1, 1186, 1243.8, 1302.4, 1362.2, 1422.9, 1484, 1546.4, 1609.6, 1673.6, 1738.4, 1804.2, 1870.7, 1937.7, 2005.4, 2143.6, 2213.4, 2284.2, 2355.3, 2427.3, 2499.7, 2572.5, 2646.2, 2720.3, 2794.7, 2870.1, 2945.4, 3021.6, 3097.7, 3174.3, 3251.7, 3329.2, 3485.3]
        self.__psi_list = [9.5266, 9.7348, 9.9243, 10.09, 10.2296, 10.3761, 10.5061, 10.618, 10.7265, 10.8285, 10.9246, 11.0156, 11.1021, 11.1846, 11.2635, 11.3392, 11.412, 11.4821, 11.5498, 11.6152, 11.6766, 11.74, 11.7995, 11.8574, 11.9137, 11.9685, 12.0219, 12.0739, 12.1246, 12.1742, 12.2226, 12.2689, 12.3161, 12.3613, 12.4056, 12.449, 12.4914, 12.5331, 12.5739, 12.6139, 12.6532, 12.7269, 12.7668, 12.8034, 12.8393, 12.8746, 12.9093, 12.9435, 12.9771, 13.0102, 13.0427, 13.0748, 13.1064, 13.1375, 13.1681, 13.1983, 13.2281, 13.2575, 13.315]
        self.__theta_psi_list = [-0.61, -0.488, -0.361, -0.231, -0.101, 0.026, 0.159, 0.285, 0.411, 0.533, 0.653, 0.769, 0.882, 0.993, 1.102, 1.208, 1.312, 1.416, 1.514, 1.612, 1.708, 1.803, 1.896, 1.988, 2.079, 2.168, 2.256, 2.342, 2.428, 2.512, 2.596, 2.678, 2.769, 2.84, 2.92, 2.999, 3.077, 3.154, 3.23, 3.305, 3.38, 3.526, 3.598, 3.67, 3.74, 3.81, 3.879, 3.947, 4.015, 4.081, 4.147, 4.213, 4.277, 4.341, 4.404, 4.467, 4.529, 4.59, 4.711]
        
        self.__T_H_inter = interpolate.interp1d(self.__temperature_list, self.__H_list)
        self.__theta_h_inter = interpolate.interp1d(self.__temperature_list, self.__theta_h_list)
        self.__psi_inter = interpolate.interp1d(self.__temperature_list, self.__psi_list)
        self.__theta_psi_inter = interpolate.interp1d(self.__temperature_list, self.__theta_psi_list)
        
        self.__psi_T_inter = interpolate.interp1d(self.__psi_list, self.__temperature_list)
        self.__H_T_inter = interpolate.interp1d(self.__H_list, self.__temperature_list)
        

    


    
    def Tdf2CP(self, T,d=0,f=0):
        """
        给定 油气比f 大气湿度d 温度T(K) 计算 定压比热容C_p
        范作民, 傅巽权. 热力过程计算与燃气表[M]. 国防工业出版社, 1987.
        """
    
        tau = T * 1e-3
        C_pa = 1e-3*sum([i * self.a_air_i[i] * (tau ** (i-1)) for i in range(9)])
        #TODO:计算湿度与油气比
        C_pa /= 1000
        return C_pa
        
    def Tdf2H(self, T,d=0,f=0):
        h_a = self.__T_H_inter(T)
        return float(h_a)

    def Tdf2H_old(self, T,d=0,f=0):
        tau = T * 1e-3
        h_a = sum([self.a_air_i[i] * (tau ** (i)) for i in range(9)]) + 301857 #偏移到参考书的零点数值
        h_a /= 1000
        return h_a

    def Tdf2theta_h(self, T,d=0,f=0):
        """
        输入温度T(K)求theta
        
        tau = T * 1e-3
        theta_h = 0
        for i in range(9):
            theta_h += self.a_theta_i[i] * (tau ** i)
        """
        #theta_h = 4e-11* (T**4) - 3e-07 * (T**3) + 0.0012 * (T**2) + 0.2589 * T - 38.681
        theta_h = self.__theta_h_inter(T)
        return float(theta_h)

    def Tdf2psi(self, T,d=0,f=0):
        """
        输入温度T(K)求psi
        tau = T * 1e-3
        psi = 1e-3 * self.a_air_i[1] * np.log(tau) + 1e-3 * sum([(i/(i-1)) * self.a_air_i[i] * (tau ** (i-1)) for i in range(2,9)]) + self.a_air_i[-1]
        """
        psi = self.__psi_inter(T)
        return float(psi)

    def Tdf2theta_psi(self, T,d=0,f=0):
        """
        输入温度T(K)求theta_psi
        """
        theta_psi = self.__theta_psi_inter(T)
        return float(theta_psi)

    def psi2T(self, psi):
        T = self.__psi_T_inter(psi)
        return float(T)

    def H2T(self, H):
        T = self.__H_T_inter(H)
        return float(T)
    '''
    def fdT2S(f=0, d, T):
        pass


    def fdT2k(f=0, d, T):
        pass


    def fdT2R(f=0, d, T):
        pass
    '''

if __name__ == "__main__":

    table = GasTable()
    temp = np.linspace(1000,2000,50)
    cp_list = []
    h_list = []
    hn_list = []
    theta_h_list = []
    psi_list = []
    theta_psi_list = []
    print("Temp \tC_p \tH \tθ_H \tψ \tθ_ψ ")
    for each in temp:
        theta_h = table.Tdf2theta_h(each)
        cp = table.Tdf2CP(each)
        h = table.Tdf2H_old(each)
        h_n = table.Tdf2H(each)
        psi = table.Tdf2psi(each)
        theta_psi = table.Tdf2theta_psi(each)
        
        theta_h_list.append(theta_h)
        cp_list.append(cp)
        h_list.append(h)
        hn_list.append(h_n)
        psi_list.append(psi)
        theta_psi_list.append(theta_psi)

        print("%f \t%f \t%f \t%f \t%f \t%f " % (each, cp, h, theta_h, psi, theta_psi))

    plt.subplot(151)
    plt.plot(temp,theta_h_list)
    plt.subplot(152)
    plt.plot(temp,cp_list)
    plt.subplot(153)
    plt.plot(temp,h_list)
    plt.plot(temp,hn_list)
    plt.subplot(154)
    plt.plot(temp,psi_list)
    plt.subplot(155)
    plt.plot(temp,theta_psi_list)
    
    plt.savefig("gas_proprity.png")
    plt.show()
