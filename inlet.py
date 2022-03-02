"""
2021-11-19
进气道参数计算
"""
import copy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import gas


class Inlet():
    def __init__(self) -> None:
        self.sigma_i = 1
        self.mach_number = 0
        

    def run(self, flow_in):
        """
        由进气总温总压 马赫数计算进气参数
        return:P_out, T_out
        """
        self.flow_in = flow_in
        self.flow_out = copy.copy(flow_in)

        #计算总压恢复系数
        self.mach_number = flow_in.c/flow_in.a
        self.__calculate_sigma_i()
        #计算出口参数
        self.flow_out.Pt = self.sigma_i*flow_in.Pt
        self.flow_out.Tt = flow_in.Tt
        self.Pi_i = self.flow_out.Pt/self.flow_in.Pt
        return self.flow_out

    def __calculate_sigma_i(self):
        """
        计算进气道总压恢复系数
        """
        if (self.mach_number <= 1.0):
            self.sigma_i = 0.97
        else:
            self.sigma_i = 0.97*(1.0-0.075*((self.mach_number-1)**1.35))
        return self.sigma_i

    def plot(self):
        point_count = 100
        mach_range = list(np.linspace(0,2,point_count))
        sigma = []
        for each in mach_range:
            self.mach_number = each
            sigma.append(self.__calculate_sigma_i())
        self.mach_number = 0
        plt.plot(mach_range, sigma)
        plt.show()

if __name__ == "__main__":
    my_inlet = Inlet()
    my_inlet.plot()