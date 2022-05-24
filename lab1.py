#####################################
###################################
#author: licabbage 
#2022.5.24
##############################
#############################

import math
import numpy as np
import matplotlib.pyplot as plt    
import pandas as pd

#####下面为最佳的一条记录

# p:min,max,setp : 3.5	3.7	0.001
# q:min,max,setp : 1.5	1.7	0.001
# r:min,max,setp : 4.1	4.3	0.001
# min_max_delta : 3.1182930595254183e-08 /3.12089538902857e-08

# p, q, r is 3.579	1.658	4.236
#二进制为
#0.0000 0000 0000 0000 0000 0000 10000110000010101001101110010001110001100100110101111
#可到达2的-25次方


class Solution:

    def __init__(self):
        #初始值，随便设的
        #最佳的一条记录如下
        # p, q, r is 3.579	1.658	4.236
        #将其赋值，运行程序直接显示
        self.optimized_p = np.float32(3.579)
        self.optimized_q = np.float32(1.658)
        self.optimized_r = np.float32(4.236)

        self.min_max_delta = np.float32(20)

        self.p_min = np.float32(3.5)
        self.p_max = np.float32(3.7)
        self.p_step = np.float32(0.01)

        self.q_min = np.float32(1.5)
        self.q_max = np.float32(1.8)
        self.q_step = np.float32(0.01)

        self.r_min = np.float32(4.1)
        self.r_max = np.float32(4.4)
        self.r_step = np.float32(0.01)

        

    def compute_min_max_delta(self):       
        min__max_delta = 20 #先给一个较大的值
        for p in np.arange(self.p_min,self.p_max,self.p_step):
            for q in np.arange(self.q_min, self.q_max,self.q_step):
                for r in np.arange(self.r_min,self.r_max,self.r_step):
                    #计算在[1:2]之间的最大误差
                    max_delta = 0.0
                    for c in np.arange(1.0,2.001,0.001,dtype=np.float32):
                        x_0 = (p*c + q) / (c+r)
                        delta =  ((x_0 - math.sqrt(c)) / (x_0 + math.sqrt(c)))**2 
                        if delta > max_delta:
                            max_delta = delta

                    if max_delta < min__max_delta:
                        print("min__max_delta update to " + str(max_delta))
                        print("current p, q, r is %.3f %.3f %.3f"%(p,q,r))
                        #将优化好的p,q,r赋值
                        self.optimized_p = p
                        self.optimized_q = q
                        self.optimized_r = r
                        min__max_delta = max_delta

        self.min_max_delta = min__max_delta
        return min__max_delta

    #画图
    def plot_min_max_delta(self):
        X_axis = np.arange(1.0,2.001,0.001,dtype=np.float32)
        #print(X_axis)
        
        x_0s = (self.optimized_p *X_axis + self.optimized_q) / (X_axis + self.optimized_r)
        #print(x_0s)

        Y_axis = ((x_0s - X_axis**0.5) / (x_0s + X_axis**0.5))**2 
        #print(Y_axis)

        plt.plot(X_axis,Y_axis)

        plt.show()

    #返回无穷范数
    def inf_norm(self,p,q,r):
        max_delta = 0.0
        for c in np.arange(1.0,2.001,0.001,dtype=np.float32):
            x_0 = (p*c + q) / (c+r)
            delta =  ((x_0 - math.sqrt(c)) / (x_0 + math.sqrt(c)))**2 
            if delta > max_delta:
                max_delta = delta
        return max_delta

    #用于一次迭代完成之后记录    
    def record(self):
        with open("record3.txt",encoding= "utf-8",mode="a") as file:
            record = "\n" +\
                    "p:min,max,setp : " + str(self.p_min) +"\t" + str(self.p_max) + "\t" + str(self.p_step) + "\n" +\
                    "q:min,max,setp : " + str(self.q_min) +"\t" + str(self.q_max) + "\t" + str(self.q_step) + "\n" +\
                    "r:min,max,setp : " + str(self.r_min) +"\t" + str(self.r_max) + "\t" + str(self.r_step) + "\n" +\
                    "min_max_delta : "  + str(self.min_max_delta) + "\n" + \
                    "p, q, r is " +str(self.optimized_p) +"\t" + str(self.optimized_q) + "\t" + str(self.optimized_r) + "\n"
                    

            file.write(record)
            print("recode finished")
        return

        




if __name__ == "__main__":

    s = Solution()

    #给定一定的参数进行迭代找结果并记录，由于已经找到，故将其注释掉不用
    # result = s.compute_min_max_delta()
    # print("min_max_result is " + str(result))
    # s.record()

    #使用当前优化后p,q,r的结果以及图示
    inf_norm = s.inf_norm(s.optimized_p,s.optimized_q,s.optimized_r)
    print("p, q, r is ",s.optimized_p,s.optimized_q,s.optimized_r)
    print("无穷范数的结果为：",inf_norm)
    print("16进制显示结果为",inf_norm.hex())
    s.plot_min_max_delta()


    
    
