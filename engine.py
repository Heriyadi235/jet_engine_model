"""
发动机部件类
2021-11-19

"""
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

import burner
import compressor
import gas
import inlet
import nozzle
import turbine


class Engine():
    def __init__(self) -> None:
        # 气体参数先定义在这里
        # TODO:从gas_property中获取气体参数
        # TODO:设计部件创建方法
        self.k_air = 1.33
        self.R_air = 287
        #构造发动机部件
        self.inlet = inlet.Inlet()
        self.comp_map = compressor.CompMap()
        self.compressor = compressor.Compressor(self.comp_map)
        self.burner = burner.Burner()
        self.turb_map = turbine.TurbineMap()
        self.turbine = turbine.Turbine(self.turb_map)
        self.nozzle = nozzle.Nozzle()

        # 定义飞行参数
        self.height = 0 #km
        self.mach_number = 0
        
        # 定义截面参数o
        # 定义初始点
        self.wf = 0.006828908418078288
        self.n = 9000 #初始转速
        self.delta_power = 0
        self.f = 0.04
        self.thrust = 0
        
        #定义初猜值
        self.X0 = [self.n, 5,2]
        #
        self.n_record = []
        self.fn_record = []
        self.Pic_record =[]
        self.Pit_record = []
        self.wa_record = []
        self.var_record = []

    def __calculate_flow(self,x_guess):
        """计算单次流道参数
        """
        self.section0 = gas.Gas()
        self.__calculate_env()

        #n_c = 8000
        #just list some guessed_value here
        #pi_c = 10
        #pi_t = 1.1
        n_c = x_guess[0]
        pi_c =x_guess[1]
        pi_t =x_guess[2]
        
        #gas params are sent to section0
        self.section2 = self.inlet.run(self.section0)
        #返回截面3 气流参数 引气参数 冷却气流参数
        self.section3, self.flow_out, self.flow_cool, self.N_c = self.compressor.run(self.section2, pi_c, n_c)
        #f = 0.03 8702 3.65 1.244
        #f = 0.04 9110 4.07 1.255
        #self.wf = self.f * self.section3.mass
        #wf = 0.006828908418078288
        self.section4 = self.burner.run(self.section3,self.wf)
        self.section5, self.N_t = self.turbine.run(self.section4, self.flow_cool, pi_t, n_c)
        self.section8 = self.nozzle.run(self.section5, self.section0)
        
        #计算推力
        self.thrust = (self.section8.c - self.section0.c) * self.section8.mass + (self.section8.P - self.section0.P) * self.nozzle.A8
        
        #计算剩余功
        self.delta_power = -self.N_t - self.N_c
        #计算残差
        e1 = (-self.N_t)/self.N_c - 1
        e2 = self.section4.mass/self.section5.mass - 1
        e3 = self.section5.mass/self.section8.mass - 1
        return [e1, e2, e3]

    def __calculate_env(self):
        """
        环境参数计算
        TODO:计算开始时从gas_property中获取气体参数
        """
        # 计算进气静压静温
        if self.height <= 11:
            #print("高度小于等于11km")
            self.section0.T = (288.15-6.5*self.height)
            self.section0.P = 101325 * ((1-self.height/44.308)**5.2553)
        else:
            #print("高度大于11km")
            self.section0.T = 216.7
            self.section0.P = 0.227*math.exp((11-self.height)/6.338)*100000

        #self.section0.P = 101325*2
        # 计算进气总压总温
        self.section0.a = math.sqrt(self.k_air*self.R_air*self.section0.T)
        self.section0.c = self.section0.a*self.mach_number
        self.section0.Pt = self.section0.P * \
            ((1+((self.k_air-1)/2)*(self.mach_number**2))
             ** (self.k_air/(self.k_air-1)))
        self.section0.Tt = self.section0.T*((1+(self.k_air-1)*(self.mach_number**2) / 2))
        
    def cal_error_stable(self,x):
        #稳态仿真的计算残差
        return self.__calculate_flow(x)

    def cal_error_dynamics(self,x):
        #动态仿真的计算残差
        x0 = [self.n, x[0], x[1]]
        y = self.__calculate_flow(x0)
        self.delta_power = y[0] 
        return [y[1], y[2]]

    def slove_dynaminc(self):
        delta_power = 0
        total_step = 2000
        for i in range(total_step):
            if i > 400:
                #self.burner.Tt4 = 850
                #self.nozzle.A8 *= 1.001
                #self.height = 0
                #self.mach_number = 0.25
                pass
            if i > 600:
                #self.burner.Tt4 = 900
                #self.nozzle.A8 *= 1
                #self.height = 0
                #self.mach_number = 0.5
                pass
            if i > 800:
                self.burner.Tt4 = 850
                #self.nozzle.A8 *= 1
                #self.height = 0.1
                #self.mach_number = 1
                pass
            if i > 1400:
                #self.burner.Tt4 = 850
                #self.nozzle.A8 *= 1
                #self.height = 0
                #self.mach_number = 1.5
                pass

            if i > 1500:
                #self.burner.Tt4 = 790
                #self.nozzle.A8 *= 1
                #self.height = 0
                #self.mach_number = 1.5
                pass

            #self.burner.Tt4 = 1500
            #self.height = i/10000
            #self.wf = 0.006828908418078288
            self.rotor_dynamics()

            x0 = [self.X0[1], self.X0[2]]
            try:
                root = fsolve(self.cal_error_dynamics, x0)
                self.X0 = [self.n, root[0], root[1]]

            except Exception:
                i = total_step
            
            self.n_record.append(self.n)
            self.fn_record.append(self.thrust/self.section8.mass)
            self.Pic_record.append(self.compressor.pi)
            self.Pit_record.append(self.turbine.pi)
            self.wa_record.append(self.section3.mass)
            self.var_record.append(self.burner.Tt4)
            pass
      
    def rotor_dynamics(self):
        #梯形法解转子动力学方程
        J_rotor = 0.3e-5 #瞎给的一个转动惯量
        tau = ((30/math.pi)**2) * (1/J_rotor)
        dt = 0.05
        n1 = self.n + (tau / self.n * self.delta_power)*dt
        self.n = self.n + 0.5* dt *((tau / self.n * self.delta_power) + (tau / n1 * self.delta_power))
        
        #print(self.delta_power)

        
#环境计算测试
if __name__ == "__main__":

    
    my_engine = Engine()
    my_engine.height = 0
    my_engine.mach_number = 0
    
    my_engine.slove_dynaminc()

    plt.subplot(3,2,1)
    plt.plot(my_engine.n_record[300:],color = 'green')
    plt.xlabel("step")
    plt.ylabel("n")
    #plt.ylim([8000,10000])
    plt.grid()
    #plt.title("shaft speed")

    plt.subplot(3,2,2)
    plt.plot(my_engine.fn_record[300:],color = 'green')
    plt.xlabel("step")
    plt.ylabel(r"$F_s$")
    #plt.ylim([0,500])
    plt.grid()
    plt.legend()

    plt.subplot(3,2,3)
    plt.plot(my_engine.Pic_record[300:],color = 'green')
    plt.xlabel("step")
    plt.ylabel(r"$\pi_c$")
    #plt.ylim([1,10])
    plt.grid()
    plt.legend()

    plt.subplot(3,2,4)
    plt.plot(my_engine.Pit_record[300:],color = 'green')
    plt.xlabel("step")
    plt.ylabel(r"$\pi_t$")
    #plt.ylim([1,4])
    plt.grid()
    plt.legend()
    
    plt.subplot(3,2,5)
    plt.plot(my_engine.wa_record[300:],color = 'green')
    plt.xlabel("step")
    plt.ylabel(r"$mass$")
    #plt.ylim([50,130])
    plt.grid()
    plt.legend()

    plt.subplot(3,2,6)
    plt.plot(my_engine.var_record[300:],color = 'green')
    plt.xlabel("step")
    plt.ylabel(r"$T_{t4}$")
    plt.grid()
    plt.legend()

    plt.show()
    
    #my_engine.compressor.map.plot()
    #my_engine.turbine.map.plot()
    
    """
    x0 = [8000, 5,1]
    my_engine.wf = 0.002
    root = fsolve(my_engine.cal_error_stable,x0)
    my_engine.X0 = root
    print(root)
    print(my_engine.cal_error_stable(root))
    """
    
