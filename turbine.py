"""
涡轮类
2021-11-19

"""
import copy
from numpy.ma.core import set_fill_value
from scipy import interpolate
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import gas
import math

from gas_table import GasTable

class TurbineMap():
    def __init__(self) -> None:
        """
        TODO:涡轮特性图如何编辑
        """
        #self.pi = [0.05, 0.1, 0.2, 0.4, 0.7, 1, 1.5, 2.5, 5, 15] # x
        self.ref_pressor = 607950
        self.ref_temperature = 1400
        self.rotating_speed = [6231, 7120, 8010, 8900, 9110, 9410, 9700, 10000,10300] # y
        #据说interpolate里有坑，数据要从小到大排
        self.pi_list = [1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.837, 3, 3.2]

        self.mass_flow = [[70.72 , 74.40 , 77.08 , 78.17 , 78.28 , 78.31 , 78.34 , 78.35 , 78.35 , 78.35 ],
                            [70.36 , 73.91 , 76.67 , 77.97 , 78.26 , 78.30 , 78.33 , 78.34 , 78.34 , 78.34],
                            [70.07 , 73.53 , 76.23 , 77.68 , 78.21 , 78.29 , 78.31 , 78.32 , 78.33 , 78.33 ],
                            [69.67 , 73.13 , 75.75 , 77.39 , 78.11 , 78.27 , 78.30 , 78.32 , 78.33 , 78.33 ],
                            [69.66 , 73.06 , 75.67 , 77.37 , 78.04 , 78.29 , 78.31 , 78.31 , 78.31 , 78.31 ],
                            [69.50 , 72.92 , 75.53 , 77.15 , 77.98 , 78.26 , 78.31 , 78.31 , 78.31 , 78.31],
                            [69.35 , 72.77 , 75.36 , 77.03 , 77.87 , 78.21 , 78.25 , 78.29 , 78.29 , 78.29 ],
                            [68.96 , 72.58 , 75.15 , 76.84 , 77.78 , 78.12 , 78.23 , 78.26 , 78.26, 78.26],
                            [68.96 , 72.58 , 75.15 , 76.84 , 77.78 , 78.12 , 78.23 , 78.26 , 78.26, 78.26]
                            ]

        self.effiency = [[0.892 , 0.870 , 0.847 , 0.825 , 0.802 , 0.779 , 0.754 , 0.726 , 0.706 , 0.682 ],
                        [0.903 , 0.896 , 0.880 , 0.866 , 0.850 , 0.832 , 0.813 , 0.790 , 0.770 , 0.744 ],
                        [0.901 , 0.903 , 0.902 , 0.894 , 0.885 , 0.873 , 0.859 , 0.840 , 0.825 , 0.808 ],
                        [0.874 , 0.892 , 0.903 , 0.907 , 0.906 , 0.901 , 0.892 , 0.880 , 0.867 , 0.851 ],
                        [0.870 , 0.890 , 0.902 , 0.908 , 0.909 , 0.905 , 0.897 , 0.887 , 0.875 , 0.859 ],
                        [0.859 , 0.884 , 0.899 , 0.907 , 0.910 , 0.909 , 0.903 , 0.893 , 0.883 , 0.870 ],
                        [0.840 , 0.871 , 0.893 , 0.904 , 0.910 , 0.910 , 0.905 , 0.897 , 0.888 , 0.875 ],
                        [0.880 , 0.895 , 0.906 , 0.912 , 0.915 , 0.915 , 0.910 , 0.900 , 0.888 , 0.872 ],
                        [0.858 , 0.886 , 0.904 , 0.912 , 0.915 , 0.914 , 0.908 , 0.898 , 0.888 , 0.875 ]
                        ]

                
        self.speed_pi2mass = interpolate.interp2d(self.pi_list, self.rotating_speed, self.mass_flow)
        self.speed_pi2eff = interpolate.interp2d(self.pi_list, self.rotating_speed, self.effiency)
        
    def get_mass(self, pi, n):
        """
        给定压比与转速，计算换算流量
        """
        correct_flow_mass = float(self.speed_pi2mass(pi, n))  #平移特性图

        #print("%f %f %f" %(pi, n, correct_flow_mass))
        return correct_flow_mass

    def get_effiency(self, pi, n):
        """
        给定压比与转速，计算效率
        """
        
        effiency = float(self.speed_pi2eff(pi, n))
        return effiency

    def plot(self):
        """
        绘制特性线
        """
        pi_min = 1.4
        line_count = 50 #插值后绘制100条
        #绘制原始流量线
        plt.subplot(221)
        for idx, n in enumerate(self.rotating_speed):
            mass = self.mass_flow[idx][:]
            plt.scatter(self.pi_list, mass, s=10)
        
        #绘制插值后流量线
        plt.subplot(222)
        for n in list(np.linspace(min(self.rotating_speed),max(self.rotating_speed),line_count)):
            pi_max = max(self.pi_list)
            mass = [] 
            pi_list_inter = list(np.linspace(pi_min,pi_max,line_count))
            for pi in pi_list_inter:
                mass.append(float(self.get_mass(pi, n)))
            plt.plot(pi_list_inter, mass)
        
        #绘制原始效率线 
        plt.subplot(223)
        for idx, n in enumerate(self.rotating_speed):          
            effi = self.effiency[idx]
            plt.scatter(self.pi_list, effi,s=10)
        
        #绘制插值后效率线
        plt.subplot(224)
        for n in list(np.linspace(min(self.rotating_speed),max(self.rotating_speed),line_count)):
            pi_max = max(self.pi_list)
            mass = [] 
            effi = []
            pi_list_inter = list(np.linspace(pi_min,pi_max,line_count))
            for pi in pi_list_inter:
                effi.append(float(self.get_effiency(pi, n)))
            plt.plot(pi_list_inter, effi)

        plt.show()
        #绘制效率线
        

class Turbine():
    def __init__(self, map):
        """
        定义特性图数据
        输入一个涡轮特性图对象
        """
        self.map = map
        self.pi = 1
        self.table_gas = GasTable()


    def run(self, flow_in, flow_cool, pi, n):
        """
        pi用的前比后，就是大于1那种
        return: 出口总压，出口总温，进口流量，出口流量，冷却流量，功率
        
        """
        if pi<1:
            pi = 1
        if n<5000:
            n = 5000
        if n>12000:
            n = 12000
        
        self.pi = pi
        self.flow_in = flow_in
        self.flow_out = copy.copy(flow_in)
        
        n_cor = n * math.sqrt(self.map.ref_temperature/flow_in.Tt)
        #n_cor = n
        w_acor = self.map.get_mass(pi, n_cor)
        ita_acor = self.map.get_effiency(pi, n_cor)
        
        w_a = w_acor * self.flow_in.Pt/ self.map.ref_pressor * math.sqrt(self.map.ref_temperature/flow_in.Tt)
        w_a = w_acor

        ita = ita_acor
        self.flow_out.Pt = self.flow_in.Pt / pi
        
        
        self.flow_in.psi = self.table_gas.Tdf2psi(self.flow_in.Tt)
        self.flow_in.H = self.table_gas.Tdf2H(self.flow_in.Tt)

        psi_out_i = self.flow_in.psi - math.log10(pi)
        Tt_out_i = self.table_gas.psi2T(psi_out_i)
        h_out_i = self.table_gas.Tdf2H(Tt_out_i)

        self.flow_out.H = self.flow_in.H - (self.flow_in.H - h_out_i)/ita
        self.flow_out.Tt = self.table_gas.H2T(self.flow_out.H)
        
        N_t = w_a * (self.flow_out.H - self.flow_in.H) 

        self.flow_out.mass = w_a + flow_cool.mass
        return self.flow_out, N_t


if __name__ == "__main__":
    #测试
    map = TurbineMap()
    map.plot()
    