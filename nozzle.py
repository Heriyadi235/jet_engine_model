"""
喷管类
2021-11-19
"""

import gas
import copy
import math

class Nozzle():
    #现在这里定义个简单收敛尾喷管
    def __init__(self) -> None:
        self.flow_in = gas.Gas()
        self.flow_out = gas.Gas()
        self.flow_env = gas.Gas()
        self.A8 = 0.24 #喷口临界面积
        self.A8 = 0.239
        self.k8 = 1.3
        self.c8 = 0
        self.Ma8 = 0
        self.R8 = 287
        self.fei = 0.97
    
    def run(self, flow_in, flow_env):
        """
        由进口总压 总温 临界面积 出口面积 背压 
        计算出口气流速度、出口静压、流量
        return:P_out, 
        """
        self.flow_env = flow_env
        self.flow_in = flow_in
        self.flow_out = copy.copy(flow_in)
        
        pi_us = self.flow_in.Pt/ self.flow_env.P #可用压比
        pi_cr = ((self.k8+1)/2) ** (self.k8/(self.k8-1))
        
        if (pi_us < pi_cr):
            #亚临界状态
            self.flow_out.P = self.flow_env.P
            self.Ma8 = math.sqrt((2 / (self.k8 - 1)) * ((self.flow_out.Pt / self.flow_out.P) ** 0.2308 - 1))
            self.c8 = self.fei * math.sqrt(2*self.flow_env.Cp*self.flow_env.Tt*(1-(self.flow_env.P/self.flow_in.Pt)**((self.flow_in.k-1)/self.flow_in.k)))
            self.flow_out.c = self.c8

        else:
            self.flow_out.P = self.flow_in.Pt / pi_cr
            self.Ma8 = 1
            self.a8 = math.sqrt(self.flow_in.k * self.flow_in.Rg * self.flow_env.Tt)
            self.c8 = self.fei * self.a8 * self.Ma8
            self.flow_out.c = self.c8
            

        flow_K = self.__calculate_K(self.flow_in.k, self.flow_in.Rg)
        self.flow_out.mass = self.__calculate_q(self.flow_in.k, self.Ma8) * flow_K * self.A8 * self.flow_in.Pt / math.sqrt(self.flow_in.Tt)
        return self.flow_out 

    def __calculate_q(self,k,Ma):
        #计算密流函数q
        lambda_gas = math.sqrt((k+1)/2*Ma*Ma/(1+(k-1)/2*Ma*Ma))
        q = ((k+1)/2)**(1/(k-1))*lambda_gas*(1-(k-1)/(k+1)*lambda_gas*lambda_gas)**(1/(k-1))
        return q   

    def __calculate_K(self,k,Rg):
        #计算流量系数K
        return math.sqrt(k/Rg*(2/(k+1))**((k+1)/(k-1)))


