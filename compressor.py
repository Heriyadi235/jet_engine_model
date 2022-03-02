"""
压气机类
2021-11-19

"""
import copy
import math
from scipy import interpolate
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import gas_table

class CompMap():
    def __init__(self) -> None:
        """
        TODO:压气机特性图如何编辑
        """
        #self.pi = [0.05, 0.1, 0.2, 0.4, 0.7, 1, 1.5, 2.5, 5, 15] # x
        self.ref_pressor = 101325
        self.ref_temperature = 288.15
        self.rotating_speed = [6000, 7000, 8000, 8500 ,9000, 9550, 10000] # y
        self.rotating_speed = [4000,5000,6000,7000,8000,9000,10000] 
        self.pi = [[1, 1.6105 , 1.8930 , 2.0221 , 2.1300 ],
                    [1, 1.7932 , 2.3439 , 2.5448 , 2.7290 ],
                    [1, 2.4765 , 3.0860 , 3.5107 , 3.6516 ],
                    [1, 3.2234 , 4.1604 , 4.4923 , 4.5149 ],
                    [1, 3.9041 , 4.8723 , 5.3049 , 5.3477 ],
                    [1, 4.4326 , 5.6781 , 6.1242 , 6.2870 ],
                    [1, 6.0456 , 6.6860 , 7.1592 , 7.2674 ]
                    ]

        self.mass_flow = [[44.033 , 41.889 , 39.797 , 37.670 , 34.667 ],
                    [51.833 , 51.597 , 50.284 , 48.040 , 42.440 ],
                    [60.163 , 58.766 , 57.902 , 53.300 , 48.453 ],
                    [70.167 , 68.069 , 66.817 , 61.984 , 57.278 ],
                    [87.000 , 85.059 , 82.649 , 78.105 , 73.200 ],
                    [94.300 , 93.500 , 92.800 , 89.932 , 84.758 ],
                    [101.000 , 100.000 , 99.733 , 98.509 , 97.128 ]
                    ] # z1
        
        self.effiency = [[0.550 , 0.655 , 0.759 , 0.757 , 0.733 ],
                        [0.593 , 0.621 , 0.794 , 0.790 , 0.740 ],
                        [0.793 , 0.826 , 0.840 , 0.793 , 0.750 ],
                        [0.800 , 0.834 , 0.852 , 0.798 , 0.750 ],
                        [0.810 , 0.844 , 0.868 , 0.821 , 0.767 ],
                        [0.790 , 0.836 , 0.863 , 0.836 , 0.785 ],
                        [0.790 , 0.844 , 0.841 , 0.818 , 0.800 ]
                        ] # z2

        
        #特性数据处理
        self.rotating_speed = np.array(self.rotating_speed)
        self.pi = np.array(self.pi)
        self.mass_flow = np.array(self.mass_flow)
        self.effiency = np.array(self.effiency)
        self.index = np.arange(len(self.pi[0])) #在等转速线上做斜线取点时的线号

        self.n_idx2pi = interpolate.interp2d(self.index, self.rotating_speed, self.pi)
        self.n_idx2wa = interpolate.interp2d(self.index, self.rotating_speed, self.mass_flow)
        self.n_idx2eff = interpolate.interp2d(self.index, self.rotating_speed, self.effiency)
        

   
    def get_mass(self, pi, n):
        """
        给定压比与转速，计算换算流量
        TODO:判断是否喘振
        """
        pi_line = [float(self.n_idx2pi(i, n)) for i in self.index] #插值得到中间压比线
        idx_line = interpolate.interp1d(pi_line, self.index) #找到输入压比所在位置
        idx = idx_line(pi) 
        flow_mass = float(self.n_idx2wa(idx, n)) #由位置查得流量
        return flow_mass

    def get_effiency(self, pi, n):
        """
        给定压比与转速，计算效率
        """
        pi_line = [float(self.n_idx2pi(i, n)) for i in self.index]
        idx_line = interpolate.interp1d(pi_line, self.index)
        idx = idx_line(pi)
        effiency = float(self.n_idx2eff(idx, n))
        return effiency

    def plot(self):
        """
        绘制特性线
        """
        line_count = 50 #插值后绘制100条
        #绘制原始流量线
        plt.subplot(221)
        for idx, n in enumerate(self.rotating_speed):
            pi_list = self.pi[idx]
            mass = self.mass_flow[idx]
            plt.scatter(mass, pi_list,s=10)
        
        #绘制插值后流量线
        
        plt.subplot(222)
        for n in list(np.linspace(min(self.rotating_speed),max(self.rotating_speed),line_count)):
            pi_max =  self.n_idx2pi(max(self.index), n)
            mass = [] 
            pi_list = list(np.linspace(1,pi_max,line_count))
            for pi in pi_list:
                mass.append(float(self.get_mass(pi, n)))
            plt.plot(mass, pi_list)
        
        
        #绘制原始效率线
        plt.subplot(223)
        for idx, n in enumerate(self.rotating_speed):
            #pi_list = 0
            mass = self.mass_flow[idx] #切个片，不画喘振那段了
            effi = self.effiency[idx]
            plt.scatter(mass, effi,s=10)
        
        #绘制插值后效率线
        plt.subplot(224)
        for n in list(np.linspace(min(self.rotating_speed),max(self.rotating_speed),line_count)):
            pi_max = self.n_idx2pi(max(self.index), n) #不喘振的最大压比
            mass = [] 
            effi = []
            pi_list = list(np.linspace(1,pi_max,line_count))
            for pi in pi_list:
                mass.append(float(self.get_mass(pi, n)))
                effi.append(float(self.get_effiency(pi, n)))
            plt.plot(mass, effi)
        
        plt.show()
        #绘制效率线
        

class Compressor():
    def __init__(self, map) -> None:
        """
        定义特性图数据
        输入一个压气机特性图对象
        """
        self.map = map
        self.table_gas = gas_table.GasTable()
        self.pi = 1

    def run(self, flow_in, pi, n):
        """
        return: 出口总压，出口总温，进口流量，出口流量，冷却流量，功率
        
        """
        if pi<1:
            pi = 1
        if pi>10:
            pi = 8
        if n<5000:
            n = 5000
        if n>12000:
            n = 12000
        
        self.pi = pi
        self.flow_in = flow_in
        self.flow_out = copy.copy(flow_in)
        n_cor = n * math.sqrt(self.map.ref_temperature/flow_in.Tt)
        w_acor = self.map.get_mass(pi, n_cor) 
        ita_acor = self.map.get_effiency(pi, n_cor)
        
        #进口物理流量，实际效率
        w_a = w_acor * self.flow_in.Pt/ self.map.ref_pressor * math.sqrt(self.map.ref_temperature/flow_in.Tt)
        ita = ita_acor

        self.flow_out.Pt = self.flow_in.Pt * pi

        self.flow_in.H = self.table_gas.Tdf2H(self.flow_in.Tt)
        self.flow_in.psi = self.table_gas.Tdf2psi(self.flow_in.Tt)
        
        psi_out_i = self.flow_in.psi+math.log10(pi)
        Tt_out_i = self.table_gas.psi2T(psi_out_i)
        h_out_i = self.table_gas.Tdf2H(Tt_out_i)
      
    
        self.flow_out.H = self.flow_in.H + (h_out_i-self.flow_in.H)/ita
        self.flow_out.Tt = self.table_gas.H2T(self.flow_out.H)
        
        self.N_c = (self.flow_out.H - self.flow_in.H) * w_a

        #分气流
        self.flow_in.mass = w_a
        self.flow_out.mass = self.flow_in.mass
        self.flow_out_beta = copy.copy(self.flow_out)
        self.flow_out_beta.mass = 0
        self.flow_cool = copy.copy(self.flow_out)
        self.flow_cool.mass = 0

        return self.flow_out, self.flow_out_beta, self.flow_cool ,self.N_c


if __name__ == "__main__":
    #测试
    map = CompMap()
    map.plot()
    