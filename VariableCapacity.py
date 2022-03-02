import math
from re import S
import gas_table

class VariableCapacity:
    def __init__(self):
        self.kair = 1.4
        self.Rair = 287
        self.Cpair = 1.005
        self.Cpgas = 1.244
        self.Pich = 4.5
        self.Itaf = 0.868
        self.Itab = 0.98
        self.Itach = 0.878
        self.Itath = 0.89
        self.Itatl = 0.91
        self.Itamh = 0.98
        self.Itaml = 0.98
        self.Itabab = 0.97
        self.sigmab = 0.97
        self.sigmabab_on = 0.96
        self.sigmabab_off = 0.98
        self.sigmam = 0.97
        self.sigma2 = 0.98
        self.sigmacr = 0.98
        self.kg = 1.3
        self.beta = 0.01
        self.delta1 = 0.05
        self.delta2 = 0.05
        self.Hu = 42900
        self.Tt7ab = 2000
        self.Ct0 =3
        self.Itamp = 0.98

        self.height = 0
        self.mach = 0
        self.B = 0.4
        self.Tt4 = 1800
        self.Pifan = 3.5
        self.sigmai = 0.97
        self.T0 = 0
        self.P0 = 0
        self.a0 = 0
        self.c0 = 0
        self.Pt0 = 0
        self.Tt0 = 0
        self.Pt2 = 0
        self.Tt2 = 0
        self.Pii = 0
        self.Pt22 = 0
        self.Tt22 = 0
        self.Lcl = 0
        self.Pt3 = 0
        self.Tt3 = 0
        self.Lch = 0
        self.f = 0
        self.Pt4 = 0
        self.Tt4a = 0
        self.Pt4a = 0
        self.Tt45 = 0
        self.Pith = 0
        self.Pitl = 0
        self.Pt45 = 0
        self.Taom1 = 0
        self.Taom2 = 0
        self.Tt45ovTt4a = 0
        self.Tt5ovTt4c = 0
        self.Tt4c = 0
        self.Pt4c = 0
        self.Tt4 = 0
        self.Pt4 = 0
        self.Tt5 = 0
        self.Pith = 0
        self.Pt5 = 0
        self.Bm = 0
        self.Cp6 = 0
        self.Tt6ovTt5 = 0
        self.Tt6 = 0
        self.Pm = 0
        self.Pt6 = 0
        self.Cp7 = 0
        self.fab = 0
        self.f0ab = 0
        self.f0 = 0
        self.Pt7ab = 0

        self.Pt7 = 0
        self.Tt7 = 0
        self.use_lavaor_nozzle = 0
        self.P9 = 0
        self.Tt9ab = 0
        self.Pt9ab = 0
        self.Ma9ab = 0
        self.T9ab = 0
        self.a9ab = 0
        self.c9ab = 0
        self.Fsab = 0
        self.sfcab = 0
        self.Tt9 = 0
        self.Pt9 = 0
        self.Ma9 = 0
        self.T9 = 0
        self.a9 = 0
        self.c9 = 0
        self.Fs = 0
        self.sfc = 0
        self.feie = 0.97
        self.table_gas = gas_table.GasTable()

    def load_design_point(self,height,mach,B,Tt4,Pifan):
        self.height = height
        self.mach = mach
        self.B = B
        self.Tt4 = Tt4
        self.Pifan = Pifan
        if (mach <= 1.0):
            self.sigmai = 0.97
        else:
            self.sigmai = 0.97 * (1.0 - 0.075 * ((mach - 1) ** 1.35))
    
    def calaulate_section_0(self):
        #环境计算
        if self.height < 11:
            self.T0 = (288.15 - 6.5 * self.height)
            self.P0 = 101325 * ((1 - self.height / 44.308) ** 5.2553)
        else:
            self.T0 = 216.7
            self.P0 = 0.227 * math.exp((11 - self.height) / 6.338) * 100000
        self.a0 = math.sqrt(self.kair * self.Rair * self.T0)
        self.c0 = self.a0 * self.mach
        self.Pt0 = self.P0 * ((1 + ((self.kair - 1) / 2) * (self.mach ** 2)) ** (self.kair / (self.kair - 1)))
        self.Tt0 = self.T0 * ((1 + (self.kair - 1) * (self.mach ** 2) / 2))

    def calaulate_section_1(self):
        #进气道计算
        self.Pt2 = self.sigmai * self.Pt0
        self.Tt2 = self.Tt0
        self.Pii = self.Pt2 / self.Pt0
   
    def calaulate_section_2(self):
        #风扇计算
        lookup_table_history = [] #返回查表的历史纪录，字符串形式
        self.Pt22 = self.Pt2 * self.Pifan
        
        self.ht2 = self.table_gas.Tdf2H(self.Tt2)
        #self.psi_t2 = 7.80201 + 0.0104143 * self.Tt2 - 0.871623e-5 * (self.Tt2 ** 2)
        self.psi_t2 = self.table_gas.Tdf2psi(self.Tt2)
        lookup_table_history.append("根据进口总温Tt2=%10.4f 查表得总焓ht2为 %10.4f " %(self.Tt2, self.ht2))
        lookup_table_history.append("根据进口总温Tt2=%10.4f 查表得熵函数psit2为 %10.4f " %(self.Tt2, self.psi_t2))
        
        self.psi_t22i = self.psi_t2+math.log10(self.Pifan) #计算理想熵函数
        self.Tt22i = self.table_gas.psi2T(self.psi_t22i)
        self.ht22i = self.table_gas.Tdf2H(self.Tt22i)
      
        lookup_table_history.append("根据理想熵函数psi_t22i=%10.4f 查表得总温度Tt22i为 %10.4f " %(self.psi_t22i, self.Tt22i))
        lookup_table_history.append("根据理想总温Tt22i=%10.4f 查表得总焓值ht22i为 %10.4f " %(self.Tt22i, self.ht22i))
        
        self.ht22 = self.ht2+(self.ht22i-self.ht2)/self.Itaf

        self.Tt22 = self.table_gas.H2T(self.ht22)
        lookup_table_history.append("根据实际风扇出口总焓ht22=%10.4f 查表得总温度Tt22为 %10.4f " %(self.ht22, self.Tt22))
        self.Lcl = self.ht22 - self.ht2
        return lookup_table_history

    def calaulate_section_3(self):
        #高压机计算
        lookup_table_history = [] #返回查表的历史纪录，字符串形式
        self.Pt3 = self.Pt22 * self.Pich
        
        #self.psi_t22 = 7.80201 + 0.0104143 * self.Tt22 - 0.871623e-5 * (self.Tt22 ** 2)
        self.psi_t22 = self.table_gas.Tdf2psi(self.Tt22)
        lookup_table_history.append("根据进口总焓值ht22=%10.4f 查表得熵函数psit22为 %10.4f " %(self.ht22, self.psi_t22))
        
        self.psi_t3i = self.psi_t22+math.log10(self.Pich)
        
        #self.Tt3i = -0.174713e4 +  0.540265e3 * self.psi_t3i - 0.495662e2 * (self.psi_t3i ** 2) + 1.60424 * (self.psi_t3i ** 3) - 0.162647 * (self.psi_t3i ** 4) + 0.0159732 * (self.psi_t3i ** 5)
        #self.ht3i = 3.95694 + 0.968169 * self.Tt3i + 0.610914e-4 * (self.Tt3i ** 2) 
        self.Tt3i = self.table_gas.psi2T(self.psi_t3i)
        self.ht3i = self.table_gas.Tdf2H(self.Tt3i)

        lookup_table_history.append("根据理想熵函数psit3i=%10.4f 查表得总温度Tt3i为 %10.4f " %(self.psi_t3i, self.Tt3i))
        lookup_table_history.append("根据理想熵函数psit3i=%10.4f 查表得总焓值ht3i为 %10.4f " %(self.psi_t3i, self.ht3i))

        self.ht3 = self.ht22+(self.ht3i-self.ht22)/self.Itach

        #self.Tt3 = -18.7787 + 1.08538 * self.ht3 - 0.108229e-3 * (self.ht3 ** 2)
        self.Tt3 = self.table_gas.H2T(self.ht3)
        lookup_table_history.append("根据实际高压机出口总焓ht=%10.4f 查表得总温度Tt3为 %10.4f " %(self.ht3, self.Tt3))
        
        self.Lch = self.ht3 - self.ht22

        return lookup_table_history

    def calaulate_section_4(self):
        #燃烧室计算
        lookup_table_history = []
        self.Pt4 = self.sigmab*self.Pt3
        self.thetaf = 0.9898
        #fa ration 5-65
        self.f = self.thetaf*(-0.0110966+0.0000197799*self.Tt4+0.00000000495727*(self.Tt4**2)+(5-0.01*self.Tt3)*(0.00258+0.0000002*self.Tt4))
        #fa ratio fix factor 5-64
        
        #self.h4a = 0.321975e3 + 0.0469811 * self.Tt4 + 0.939525e-3 * (self.Tt4**2) - 0.262836e-6 * (self.Tt4**3)
        self.h4a = self.table_gas.Tdf2H(self.Tt4)
        #self.thetah = -0.160747e3 + 0.884484 * self.Tt4 + 0.175335e-3 * (self.Tt4**2) + 0.408376e-6 * (self.Tt4**3) - 0.140527e-9 * (self.Tt4**4)
        self.thetah = self.table_gas.Tdf2theta_h(self.Tt4)
        lookup_table_history.append("根据Tt4=%10.4f 查表得h4a为 %10.4f " %(self.Tt4, self.h4a))
        lookup_table_history.append("根据Tt4=%10.4f 查表得thetah为 %10.4f " %(self.Tt4, self.thetah))
        
        self.ht4g = self.h4a+self.f/(1+self.f)*self.thetah
        return lookup_table_history

    def calaulate_section_45(self):
        #高压涡轮计算
        lookup_table_history = []
        self.section_45_theta_history = []
        self.section_45_htaa_history = []
        self.Pt4a = self.Pt4
        self.ht3a = self.ht3
        
        self.ht4ag = (self.ht4g*(1-self.beta-self.delta1-self.delta2)*(1+self.f)+self.ht3a*self.delta1)/((1-self.beta-self.delta1-self.delta2)*(1+self.f)+self.delta1)
        
        #lookup_table_history.append("根据高压转子燃气进口焓值 %10.4f 计算温度" %(self.ht4ag))
        #燃气焓值对应温度迭代
        self.ht4aa = self.ht4ag

        while 1:
            #theta = -206.027 + 0.949191 * self.ht4aa + 0.348278e-3 * (self.ht4aa**2)
            temp = self.table_gas.H2T(self.ht4aa)
            theta = self.table_gas.Tdf2theta_h(temp)
            self.section_45_theta_history.append(theta)
            lookup_table_history.append("根据ht4aa=%10.4f 查表得thetah为 %10.4f" %(self.ht4aa,theta))
            de = abs((self.ht4ag - self.f/(1+self.f)* theta)-self.ht4aa)
            self.ht4aa = self.ht4ag - self.f/(1+self.f)* theta
            self.section_45_htaa_history.append(self.ht4aa)

            if de<0.0001:
                break

        #self.Tt4aa = 27.1498 + 0.979141*self.ht4aa - 0.464182e-4 * (self.ht4aa**2)
        self.Tt4aa = self.table_gas.H2T(self.ht4aa)
        lookup_table_history.append("根据ht4aa=%10.4f 查表得Tt4aa为 %10.4f" %(self.ht4aa,self.Tt4aa))

        #高压涡轮出口总温
        self.Lth = self.Lch/(((1-self.beta-self.delta1-self.delta2)*(1+self.f)+ self.delta1)* self.Itamh)
        self.ht45g = self.ht4ag - self.Lth

        self.ht45a = self.ht45g
        while 1:
            temp = self.table_gas.H2T(self.ht45a)
            theta = self.table_gas.Tdf2theta_h(temp)

            lookup_table_history.append("根据ht45a=%10.4f 查表得theta为 %10.4f" %(self.ht45a, theta))
            
            de = abs((self.ht45g - self.f/(1+self.f)*theta)-self.ht45a)
            self.ht45a = self.ht45g - self.f/(1+self.f)*theta
            self.section_45_theta_history.append(theta)
            self.section_45_htaa_history.append(self.ht45a)

            if de<0.0001:
                break
        
        self.Tt45a = self.table_gas.H2T(self.ht45a)
        lookup_table_history.append("根据ht45a=%10.4f 查表得Tt45a为 %10.4f" %(self.ht45a,self.Tt4aa))
        self.Tt45 = self.Tt45a
        #高压涡轮膨胀比及出口压力计算
        self.ht45gi = self.ht4ag -(self.Lth)/(self.Itath)
        self.ht45ai = self.ht45gi
        while 1:
            #theta = -206.027 + 0.949191 * self.ht45ai + 0.348278e-3 * (self.ht45ai**2)
            temp = self.table_gas.H2T(self.ht45ai)
            theta = self.table_gas.Tdf2theta_h(temp)

            lookup_table_history.append("根据ht45ai=%10.4f 查表得theta为 %10.4f" %(self.ht45ai, theta))
            
            de = abs((self.ht45gi - self.f/(1+self.f)*theta)-self.ht45ai)
            self.ht45ai = self.ht45gi - self.f/(1+self.f)*theta
            self.section_45_theta_history.append(theta)
            self.section_45_htaa_history.append(self.ht45ai)
            if de<0.0001:
                break

        #self.Tt45ai = 27.1498 + 0.979141*self.ht45ai - 0.464182e-4 * (self.ht45ai**2)
        self.Tt45ai = self.table_gas.H2T(self.ht45ai)
        lookup_table_history.append("根据ht45ai=%10.4f 查表得Tt45ai为 %10.4f" %(self.ht45ai,self.Tt45ai))
        
        #计算熵函数
        self.psi4a = self.table_gas.Tdf2psi(self.Tt4)
        self.theta_psi_4a = self.table_gas.Tdf2theta_psi(self.Tt4)
        
        lookup_table_history.append("根据Tt4a=%10.4f 查表得psi4a为 %10.4f" %(self.Tt4,self.psi4a))
        lookup_table_history.append("根据Tt4a=%10.4f 查表得theta_psi4a为 %10.4f" %(self.Tt4,self.theta_psi_4a))

        self.psi45ai = self.table_gas.Tdf2psi(self.Tt45ai)
        self.theta_psi_45ai = self.table_gas.Tdf2theta_psi(self.Tt45ai)
        
        lookup_table_history.append("根据Tt45ai=%10.4f 查表得psi45ai为 %10.4f" %(self.Tt45ai,self.psi45ai))
        lookup_table_history.append("根据Tt45ai=%10.4f 查表得theta_psi45ai为 %10.4f" %(self.Tt45ai,self.theta_psi_45ai))

        self.psi4ag = self.psi4a + self.f/(1+self.f)*self.theta_psi_4a
        self.psi45gi = self.psi45ai + self.f/(1+self.f)*self.theta_psi_45ai
        
        self.Pith = 10**(self.psi4ag-self.psi45gi) #584
        self.Pt45 = self.Pt4a/self.Pith
        return lookup_table_history

    def calaulate_section_5(self):
        #低压涡轮计算
        lookup_table_history = []
        self.section_5_theta_history = []
        self.section_5_htaa_history = []
        self.Pt45a = self.Pt45
        self.ht3a = self.ht3
        
        self.ht4cg = (self.ht45g*(1-self.beta-self.delta1- self.delta2)*(1+self.f)+self.ht3a*(self.delta2))/((1-self.beta-self.delta1 - self.delta2)*(1+self.f)+(self.delta2))
        #lookup_table_history.append("根据低压转子燃气进口焓值 %10.4f 计算温度" %(self.ht4cg))
        #燃气焓值对应温度迭代
        self.ht4ca = self.ht4cg

        while 1:
            
            temp = self.table_gas.H2T(self.ht4ca)
            theta = self.table_gas.Tdf2theta_h(temp)
            lookup_table_history.append("根据ht4ca=%10.4f 查表得thetah为 %10.4f" %(self.ht4ca,theta))
            de = abs((self.ht4cg - self.f/(1+self.f)* theta)-self.ht4ca)
            self.ht4ca = self.ht4cg - self.f/(1+self.f)* theta
            
            self.section_5_theta_history.append(theta)
            self.section_5_htaa_history.append(self.ht4ca)

            if de<0.0001:
                break

        #self.Tt4ca = 27.1498 + 0.979141*self.ht4ca - 0.464182e-4 * (self.ht4ca**2)
        self.Tt4ca = self.table_gas.H2T(self.ht4ca)
        lookup_table_history.append("根据ht4ca %10.4f 查表得Tt4ca %10.4f" %(self.ht4ca,self.Tt4ca))

        #低压涡轮出口总温
        #TODO:下式是否正确（效率与流量）？
        #self.Ltl = self.Lcl *(1-self.beta-self.delta1-self.delta2)/((1-self.beta-self.delta1-self.delta2)*(1+self.f)* self.Itaml)
        self.Ltl = (self.Lcl+self.Ct0/self.Itamp)*(1+self.B)/(((1-self.beta-self.delta1-self.delta2)*(1+self.f)+ self.delta2+ self.delta1)* self.Itaml)
        self.ht5g = self.ht4cg - self.Ltl

        self.ht5a = self.ht5g

        while 1:
            #theta = -206.027 + 0.949191 * self.ht5a + 0.348278e-3 * (self.ht5a**2)
            temp = self.table_gas.H2T(self.ht5a)
            theta = self.table_gas.Tdf2theta_h(temp)
            lookup_table_history.append("根据ht5a %10.4f 查表得theta %10.4f" %(self.ht5a, theta))
            
            de = abs((self.ht5g - self.f/(1+self.f)*theta)-self.ht5a)
            self.ht5a = self.ht5g - self.f/(1+self.f)*theta

            self.section_5_theta_history.append(theta)
            self.section_5_htaa_history.append(self.ht5a)

            if de<0.0001:
                break
        
        #self.Tt5a = 27.1498 + 0.979141*self.ht5a - 0.464182e-4 * (self.ht5a**2)
        self.Tt5a = self.table_gas.H2T(self.ht5a)
        lookup_table_history.append("根据ht5a=%10.4f 查表得Tt5a为 %10.4f" %(self.ht5a,self.Tt5a))
        self.Tt5 = self.Tt5a
        #低压涡轮膨胀比及出口压力计算
        self.ht5gi = self.ht4cg -(self.Ltl)/(self.Itatl)

        #迭代ht5ai
        self.ht5ai = self.ht5gi
        while 1:
            #theta = -206.027 + 0.949191 * self.ht5ai + 0.348278e-3 * (self.ht5ai**2)
            temp = self.table_gas.H2T(self.ht5ai)
            theta = self.table_gas.Tdf2theta_h(temp)
            lookup_table_history.append("根据ht5ai=%10.4f 查表得thetah为 %10.4f" %(self.ht5ai, theta))
            
            de = abs((self.ht5gi - self.f/(1+self.f)*theta)-self.ht5ai)
            self.ht5ai = self.ht5gi - self.f/(1+self.f)*theta
            self.section_5_theta_history.append(theta)
            self.section_5_htaa_history.append(self.ht5ai)
            if de<0.0001:
                break

        self.Tt5ai = self.table_gas.H2T(self.ht5ai)
        lookup_table_history.append("根据ht5ai=%10.4f 查表得Tt5ai为 %10.4f" %(self.ht5ai,self.Tt5ai))

        #计算熵函数
        self.psi45a = self.table_gas.Tdf2psi(self.Tt45a)
        self.theta_psi_45a = self.table_gas.Tdf2theta_psi(self.Tt45a)
        lookup_table_history.append("根据Tt45a=%10.4f 查表得psi45a为 %10.4f" %(self.Tt45a,self.psi45a))
        lookup_table_history.append("根据Tt45a=%10.4f 查表得theta_psi45a为 %10.4f" %(self.Tt45a,self.theta_psi_45a))

        self.psi5ai = self.table_gas.Tdf2psi(self.Tt5ai)
        self.theta_psi_5ai = self.table_gas.Tdf2theta_psi(self.Tt5ai)

        lookup_table_history.append("根据Tt5ai=%10.4f 查表得psi5ai为 %10.4f" %(self.Tt5ai,self.psi5ai))
        lookup_table_history.append("根据Tt5ai=%10.4f 查表得theta_psi5ai为 %10.4f" %(self.Tt5ai,self.theta_psi_5ai))

        self.psi45ag = self.psi45a + self.f/(1+self.f)*self.theta_psi_45a
        self.psi5gi = self.psi5ai + self.f/(1+self.f)*self.theta_psi_5ai
        self.Pitl = 10**(self.psi45ag-self.psi5gi) #584
        self.Pt5 = self.Pt45a/self.Pitl
        
        return lookup_table_history

    def calaulate_section_6(self):
        lookup_table_history = []
        self.section_6_theta_history = []
        self.section_6_ha_history = []
        self.ht6g = (((1-self.delta1-self.delta2-self.beta)*(1+self.f)+self.delta1+self.delta2)*self.ht5g+self.B*self.ht22)/((1-self.delta1-self.delta2-self.beta)*(1+self.f)+self.delta1+self.delta2+self.B)
        self.ht6g = (((1-self.delta1-self.delta2-self.beta)*(1+self.f)+self.delta1)*self.ht5g+self.B*self.ht22)/((1-self.delta1-self.delta2-self.beta)*(1+self.f)+self.delta1+self.B)

        self.f6 = (self.f*(1-self.delta1-self.delta2-self.beta))/(1+self.B-self.beta)
        
        self.ht6a = self.ht6g
        while 1:
            temp = self.table_gas.H2T(self.ht6a)
            theta = self.table_gas.Tdf2theta_h(temp)
            lookup_table_history.append("根据ht6a=%10.4f 查表得thetah为 %10.4f" %(self.ht6a, theta))
            
            de = abs((self.ht6g - self.f6/(1+self.f6)*float(theta))-self.ht6a)
            self.ht6a = self.ht6g - self.f6/(1+self.f6)*float(theta)
            self.section_6_theta_history.append(theta)
            self.section_6_ha_history.append(self.ht6a)

            if de<0.0001:
                break
        self.Tt6 = self.table_gas.H2T(self.ht6a)
        lookup_table_history.append("根据ht6a=%10.4f 查表得Tt6为 %10.4f" %(self.ht6a, self.Tt6))

        self.Bm = self.B / ((1 - self.delta1 - self.delta2 - self.beta) * (1 + self.f) + self.delta1 + self.delta2)
        self.Pm = (self.Pt5+self.Bm*self.sigma2*self.Pt22)/(1+self.Bm)
        self.Pt6 = self.Pm*self.sigmam
        return lookup_table_history

    def calaulate_section_7(self):
        lookup_table_history = []
        self.section_7_fab_history = []
        self.section_7_ha_history = []

        #加力燃烧室计算
        #加力打开时
        self.Pt7ab = self.sigmabab_on*self.Pt6
        #Tt7 = 2000时的htab 与 theta
        self.htab = 2251.07
        self.theta2k = 3485.3
        self.section_7_ha_history.append(self.htab)
        
        self.htabg = self.htab
        
        while 1:
            lookup_table_history.append("iteration")
            self.fab = (1+(self.f*(1-self.delta1-self.delta2-self.beta))/(1+self.B-self.beta))*(self.htabg-self.ht6g)/(self.Itabab*self.Hu-self.htabg)
            de = abs(self.htab+(self.fab/(1+self.fab))*self.theta2k - self.htabg)
            self.htabg = self.htab+(self.fab/(1+self.fab))*self.theta2k
            self.section_7_fab_history.append(self.fab)
            self.section_7_ha_history.append(self.htabg)
            if de<1e-6:
                break
        
        self.f0ab = ((1-self.beta-self.delta1-self.delta2)*self.f+(1+self.B-self.beta)*self.fab)/(1+self.B)
        self.Pt7 =self.sigmabab_off*self.Pt6
        self.Tt7 = self.Tt6
        self.f0 =  ((1-self.beta-self.delta1-self.delta2)*self.f)/(1+self.B)
        return lookup_table_history

    def calaulate_section_8(self):
        #尾喷管计算 假设完全膨胀

        self.P9 = self.P0
        self.Tt9ab = self.Tt7ab
        self.Pt9ab = self.sigmacr * self.Pt7ab

        self.Tt9 = self.Tt7
        self.Pt9 = self.sigmacr * self.Pt7

        if((self.Pt9/self.P0)>5):
            self.use_lavaor_nozzle = True
            self.Ma9ab = math.sqrt((2 / (self.kg - 1)) * ((self.Pt9ab / self.P9) ** 0.2308 - 1))
            self.Ma9 = math.sqrt((2 / (self.kg - 1)) * ((self.Pt9 / self.P9) ** 0.2308 - 1))
           
        else:
            if((self.Pt9/self.P0)<1.89):
                print("P wrong !!!")
            self.use_lavaor_nozzle = False
            self.Ma9ab = 1
            self.Ma9 = 1
        

        self.T9ab = self.Tt9ab * ((1 + ((self.kg - 1) * (self.Ma9ab ** 2)) / (2)) ** (-1))
        self.a9ab = math.sqrt(self.kg * self.Rair * self.T9ab)
        self.c9aba = self.a9ab * self.Ma9ab
        self.c9ab = self.c9aba * self.feie


        self.Fsab = (1 + self.f0ab - (self.beta) / (1 + self.B)) * (self.c9ab + self.Rair * self.T9ab * (1 - self.P0 / self.P9)) - self.c0
        self.sfcab = 3600 * self.f0ab / self.Fsab

        self.T9 = self.Tt9 * ((1 + ((self.kg - 1) * (self.Ma9 ** 2)) / (2)) ** (-1))
        self.a9 = math.sqrt(self.kg * self.Rair * self.T9)
        self.c9a = self.a9 * self.Ma9
        self.c9 = self.c9a * self.feie
        self.Fs = (1 + self.f0 - (self.beta) / (1 + self.B)) * (self.c9 + self.Rair * self.T9 * (1 - self.P0 / self.P9)) - self.c0
        self.sfc = 3600 * self.f0 / self.Fs

        
