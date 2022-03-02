"""
燃烧室
2021-11-19
"""
import copy
import gas_table

class Burner():
    def __init__(self) -> None:
        self.table_gas = gas_table.GasTable()
        self.ita_b = 0.9 #燃烧效率
        self.Hu = 43124
        self.Tt4 = 800 #800~1250
        
    
    def run(self, flow_in, w_f):
        """
        根据 进口总压 总温 流量 燃油流量
        计算 出口总压 总温 燃气流量 燃烧产生的能量(?)
        """
        self.flow_in = flow_in
        self.flow_out = copy.copy(flow_in)
        self.flow_out.mass += w_f #计算总流量
        
        self.flow_out.H = (self.flow_in.mass * self.flow_in.H + w_f * self.Hu * self.ita_b)/(self.flow_out.mass)
        self.flow_out.Tt = self.flow_out.H / 1.005

        #固定一下出口温度
        self.flow_out.Tt = self.Tt4
        self.flow_out.H = 1.1 * self.flow_out.Tt
        self.flow_out.H = self.table_gas.Tdf2H(self.flow_out.Tt)

        return self.flow_out