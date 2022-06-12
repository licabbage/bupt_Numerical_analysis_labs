#####################################
###################################
#author: licabbage 
#2022.5.24
##############################
#############################

import numpy as np
import time

class Integration:
    def __init__(self) -> None:
        self.verbose = False #verbose为true可在控制台显示中间的结果

        #用于实现高斯-勒让德公式
        self.x = np.array([0.00, 
        -0.8360311073266358, 0.8360311073266358, -0.9681602395076261, 0.9681602395076261,
        -0.3242534234038089, 0.3242534234038089, -0.6133714327005904, 0.6133714327005904],
        dtype=np.float64)

        self.A = np.array([0.3302393550012598,
        0.1806481606948574, 0.1806481606948574, 0.0812743883615744, 0.0812743883615744,
        0.3123470770400029, 0.3123470770400029, 0.2606106964029354, 0.2606106964029354],
        dtype=np.float64)

    #################
    #高斯-切比雪夫 Ⅰ型
    #################
    def gauss_ch1(self, func,n:np.int32) -> np.float64:
        x =np.cos ((np.arange(1,n+1,1,dtype=np.float64)*2 -1 )  * np.pi / (2*n))

        if self.verbose:
            print("选的 %i 个点为:"%n)
            print(x)

        f_x = func(x)
        if self.verbose:
            print("选的 %i 个点对应的f为:"%n)
            print(f_x)

        result = np.sum(f_x) * np.pi / n
        if self.verbose:
            print("积分结果为：",result)

        return result

    #################
    #高斯-切比雪夫 Ⅱ型
    #################
    def gauss_ch2(self, func,n:np.int32) -> np.float64:
        points = np.arange(1,n+1,1,dtype=np.float64) * np.pi / (n+1)
        x = np.cos(points)
        if self.verbose:
            print("选的 %i 个点为:"%n)
            print(x)
        
        f_x = func(x)
        if self.verbose:
            print("选的 %i 个点对应的f为:"%n)
            print(f_x)

        sin_points = np.sin(points)
        sin_2_points = np.power(sin_points,2)
        if self.verbose:
            print("选的 %i 个点对应的sin^2为:"%n)
            print(sin_2_points)

        result = np.multiply(sin_2_points,f_x).sum() * np.pi /(n+1)
        if self.verbose:
            print("积分结果为：",result)
        return result

    #################
    #高斯-勒让德（gause-legendre) 9点,[-1,1]区间
    #################
    def gauss_leg_9(self, func) -> np.float64:
        # x = np.array([0.00, 
        # -0.8360311073266358, 0.8360311073266358, -0.9681602395076261, 0.9681602395076261,
        # -0.3242534234038089, 0.3242534234038089, -0.6133714327005904, 0.6133714327005904],
        # dtype=np.float64)

        # A = np.array([0.3302393550012598,
        # 0.1806481606948574, 0.1806481606948574, 0.0812743883615744, 0.0812743883615744,
        # 0.3123470770400029, 0.3123470770400029, 0.2606106964029354, 0.2606106964029354],
        # dtype=np.float64)

        f_x = func(self.x)

        result = np.multiply(f_x, self.A).sum()
        if self.verbose:
            print("积分结果为：",result)
        return result


    #################
    #高斯-勒让德（gause-legendre) 9点,任意区间
    #################
    def gauss_leg_9_common_range(self, func,a:np.float64, b:np.float64) -> np.float64:
        x = np.array([0.00, 
        -0.8360311073266358, 0.8360311073266358, -0.9681602395076261, 0.9681602395076261,
        -0.3242534234038089, 0.3242534234038089, -0.6133714327005904, 0.6133714327005904],
        dtype=np.float64)

        A = np.array([0.3302393550012598,
        0.1806481606948574, 0.1806481606948574, 0.0812743883615744, 0.0812743883615744,
        0.3123470770400029, 0.3123470770400029, 0.2606106964029354, 0.2606106964029354],
        dtype=np.float64)

        points = 0.5* (a+b) + 0.5*(b-a) *x
        result = 0.5*(b-a) * np.sum(A * func(points) )
        return result

    ############
    #一般区间的2点 高斯-勒让德（gause-legendre） 公式  func为函数,[a,b]为一般区间
    ############
    def gauss_leg_two_point(self,func,a:np.float64,b:np.float64)->np.float64:
        tmp1 = (b-a)/(2*np.sqrt(3))
        tmp2 = (a+b)/2
        point1 = -tmp1 + tmp2
        point2 = tmp1 + tmp2

        if self.verbose :
            print("point1 is %.3f"%point1 , " point2 is %.3f"%point2)

        result = (func(point1) + func(point2)) * (b-a) / 2
        if self.verbose :
            print("积分结果为:",result)
        return result  

    ####################
    #复化的 逐次减半 2点高斯 高斯-勒让德（gause-legendre） 公式，当前结果与上一结果误差小于一定后返回
    ###################
    def comp_gauss_leg(self,func,a:np.float64,b:np.float64,error_threshold:np.float64 = 0.000001) -> np.float64:
        pre_result = self.gauss_leg_two_point(func,a,b)

        curr_n = int(2)
        curr_result = self.gauss_leg_two_point(func,a,(a+b)/2) + self.gauss_leg_two_point(func,(a+b)/2,b)
        error = np.abs(pre_result - curr_result)

        if self.verbose:
            print("curr_n is %d "%curr_n, "error is %.5f "%error, "curr_result is %.5f "%curr_result)

        while error > error_threshold:
            pre_result = curr_result
            curr_result =0
            curr_n *= 2
            h = (b-a) /curr_n
            points = np.arange(a,b+h,h)
            assert(points.size == curr_n+1)
            for i in range(0,curr_n):
                curr_result += self.gauss_leg_two_point(func,points[i],points[i+1])

            error = np.abs(pre_result-curr_result)
            if self.verbose:
                print("curr_n is %d "%curr_n, "error is %.5f "%error, "curr_result is %.5f "%curr_result)

        return curr_result

    ############################
    #逐次减半复化梯形公式
    ############################
    def comp_trep(self,func,a:np.float64,b:np.float64, error_threshold:np.float64 = 0.000001)-> np.float64 : 
        n_divide_iter = 0
        n = 1
        h = (b-a)/n
        T_n = 0.5*h*(func(a) + func(b)) 

        if self.verbose:
            print("n_divide_iter = %i"%n_divide_iter, "T%i  = "%n, T_n)

        k = np.arange(0,n,1,dtype=np.float64)
        x_k_plus_1d2 = a + (k + 0.5) * h
        f_x_k_plus_1d2 = func(x_k_plus_1d2)
        T_2n = 0.5* T_n + 0.5*h* (f_x_k_plus_1d2.sum())

        n_divide_iter +=1
        n*=2

        if self.verbose:
            print("n_divide_iter = %i"%n_divide_iter, "T%i  = "%n, T_2n)

        error = np.abs(T_2n - T_n)
        if self.verbose:
            print("error is ", error)
        
        while error > error_threshold:
            h *= 0.5
            T_n = T_2n
            
            k = np.arange(0,n,1,dtype=np.float64)
            x_k_plus_1d2 = a + (k + 0.5) * h
            f_x_k_plus_1d2 = func(x_k_plus_1d2)

            T_2n = 0.5* T_n + 0.5*h* (f_x_k_plus_1d2.sum())
            n_divide_iter +=1
            n*=2

            if self.verbose:
                print("n_divide_iter = %i"%n_divide_iter, "T%i  = "%n, T_2n)

            error = np.abs(T_2n - T_n)
            if self.verbose:
                print("error is ", error)

        return T_2n

    ############################
    #Romberg求积公式
    ###########################
    def romberg(self,func,a:np.float64,b:np.float64,error_threshold:np.float64 = 0.000001)->np.float64:
        #一些更新用到的权重
        S_weight1 = np.float64(4.0/3.0)
        S_weight2 = -np.float64(1.0/3.0)
        C_weight1 = np.float64(16.0/15.0)
        C_weight2 = -np.float64(1.0/15.0)
        R_weight1 = np.float64(64.0/63.0)
        R_weight2 = -np.float64(1.0/63.0)
        W_weight1 = np.float64(256.0/255.0)
        W_weight2 = -np.float64(1.0/255.0)

        #初始化 S_nd2 意思为 s_n/2(下标为n/2),其他同理
        n = 1
        h = (b-a)/n
        T_n = 0.5 * (b-a) * (func(a) + func(b))
        T_2n = np.nan
        S_nd2 = np.nan
        S_n = np.nan
        C_nd4 = np.nan
        C_nd2 = np.nan
        R_nd8 = np.nan
        R_nd4 = np.nan
        W_nd16 = np.nan
        W_nd8 = np.nan

        if self.verbose:
            print("T%i is"%(n),T_n)
        #循环
        error = 20.0 #初始的一个足够大的error
        while True:
            #先更新T_2n
            k = np.arange(0,n,1,dtype=np.float64)
            x_k_plus_1d2 = a + (k + 0.5) * h
            f_x_k_plus_1d2 = func(x_k_plus_1d2)
            T_2n = 0.5* T_n + 0.5*h* (f_x_k_plus_1d2.sum())
            if self.verbose:
                print("T%i is"%(n*2),T_2n)
    
            if W_nd16 is not np.nan :
                error = np.abs(T_2n - W_nd16)
                if self.verbose:
                    print("error between T%i"%(n*2) ," and W%d is "%(n/16), "%.10f"%error)
            elif R_nd8 is not np.nan:
                error = np.abs(T_2n - R_nd8)
                if self.verbose:
                    print("error between T%i"%(n*2) ," and R%d is "%(n/8), "%.10f"%error)
            elif C_nd4 is not np.nan:
                error = np.abs(T_2n - C_nd4)
                if self.verbose:
                    print("error between T%i"%(n*2) ," and C%d is "%(n/4), "%.10f"%error)
            elif S_nd2 is not np.nan:
                error = np.abs(T_2n - S_nd2)
                if self.verbose:
                    print("error between T%i"%(n*2) ," and S%d is "%(n/2), "%.10f"%error)
            else:
                error = np.abs(T_2n - T_n)
                if self.verbose:
                    print("error between T%i"%(n*2) ," and T%d is "%n, "%.10f"%error)
            
            if(error < error_threshold):
                return T_2n

            #用 T_n 和 T_2n 更新 S_n
            if(T_n is not np.nan and T_2n is not np.nan):
                S_n = S_weight1 * T_2n + S_weight2 * T_n
                if self.verbose:
                    print("S%i is"%n,S_n)
                error = np.abs(S_n - T_2n)
                if self.verbose:
                    print("error between S%i"%(n) ," and T%d is "%(n*2), "%.10f"%error)
                if(error < error_threshold):
                    return S_n
            #用 S_nd2 和S_n 更新 C_nd2
            if(S_n is not np.nan and S_nd2 is not np.nan):
                C_nd2 = C_weight1 * S_n + C_weight2 * S_nd2
                if self.verbose:
                    print("C%i is"%(n/2),C_nd2)
                error = np.abs(C_nd2 - S_n)
                if self.verbose:
                    print("error between c%i"%(n/2) ," and S%d is "%(n), "%.10f"%error)
                if(error < error_threshold):
                    return C_nd2
                
            #用 C_nd4 和 C_nd2 更新 R_nd4
            if(C_nd4 is not np.nan and C_nd2 is not np.nan):
                R_nd4 = R_weight1 * C_nd2 + R_weight2 * C_nd4
                if self.verbose:
                    print("R%i is"%(n/4),R_nd4)
                error = np.abs(R_nd4 - C_nd2)
                if self.verbose:
                    print("error between R%i"%(n/4) ," and C%d is "%(n/2), "%.10f"%error)
                if(error < error_threshold):
                    return R_nd4

            #用 R_nd8 和 R_nd4 更新 W_nd8
            if(R_nd8 is not np.nan and R_nd4 is not np.nan):
                W_nd8 = W_weight1 * R_nd4 + W_weight2 * R_nd8
                if self.verbose:
                    print("W%i is"%(n/8),W_nd8)
                error = np.abs(W_nd8 - R_nd4)
                if self.verbose:
                    print("error between W%i"%(n/8) ," and R%d is "%(n/4), "%.10f"%error)
                if(error < error_threshold):
                    return W_nd8


            T_n = T_2n
            S_nd2 = S_n
            C_nd4 = C_nd2
            R_nd8 = R_nd4
            W_nd16 = W_nd8
            n *=2
            h *=0.5
        
        return

class normal_functions:
    #f(x) = x
    def x(input:np.ndarray):
        return input

    #f(x) = 1
    def one(input:np.ndarray):
        if type(input) == float or type(input) == int:
            return np.float64(1.0)
        return np.ones(input.shape,dtype=np.float64)

    #f(x) = x^2
    def x_power_2(input:np.ndarray):
        return np.power(input,2)

    #f(x) = sin(x) / x
    def sinx_divide_x(input:np.ndarray):
        if type(input) == float or type(input) == int:
            input = np.float64(input)
            if input == 0:
                return np.float64(1.0)
            return np.sin(input) / input
        
        result = np.sin(input) / input
        for element in result:
            if element==np.nan:
                element = 1.0
        return result

    #f(x) = x * exp(x)
    def x_mul_exp_x(input:np.ndarray):
        return np.multiply(input,np.exp(input))

    #f(x) = exp(x)
    def exp_x(input:np.ndarray):
        return np.exp(input)

    #f(x) = exp(x) * sqrt(1-x^2)
    def exp_x_mul_sqrt_1_sub_x_pow_2(input:np.ndarray):
        exp_x = np.exp(input)
        sqrt_1_sub_x_pow_2 = np.sqrt(1-np.power(input,2))
        return np.multiply(exp_x,sqrt_1_sub_x_pow_2)

    #f(x) = sin_x
    def sin_x(input:np.ndarray):
        return np.sin(input)

    #测试函数
    def test_those_methods(n_test_loops:int = 100):
        tester  = Integration()
        print("#########################\n计算exp(x) * sqrt(1-x^2),在[-1,1]区间上的积分\n#######################")

        t_begin = time.time()
        for i in range(0,n_test_loops):
            gauss_chi_1_result = tester.gauss_ch1(normal_functions.x_mul_exp_x,5)
        t_end = time.time()
        print("gauss_chi_1 result is ", gauss_chi_1_result, "\tspend time: ", (t_end-t_begin)*1000 /n_test_loops, "ms")

        t_begin = time.time()
        for i in range(0,n_test_loops):
            gauss_chi_2_result = tester.gauss_ch2(normal_functions.exp_x,5)
        t_end = time.time()
        print("gauss_chi_2 result is ", gauss_chi_2_result, "\tspend time: ", (t_end-t_begin)*1000 / n_test_loops, "ms")

        t_begin = time.time()
        for i in range(0,n_test_loops):
            gauss_leg_9_result = tester.gauss_leg_9(normal_functions.exp_x_mul_sqrt_1_sub_x_pow_2)
        t_end = time.time()
        print("gauss_leg_9 result is ", gauss_leg_9_result, "\tspend time: ", (t_end-t_begin)*1000 / n_test_loops, "ms")

        t_begin = time.time()
        for i in range(0,n_test_loops):
            comp_gauss_leg_result = tester.comp_gauss_leg(normal_functions.exp_x_mul_sqrt_1_sub_x_pow_2,-1.0,1.0)
        t_end = time.time()
        print("comp_gauss_leg result is ", comp_gauss_leg_result, "\tspend time: ", (t_end-t_begin)*1000 / n_test_loops, "ms")

        t_begin = time.time()
        for i in range(0,n_test_loops):
            comp_trep_result = tester.comp_trep(normal_functions.exp_x_mul_sqrt_1_sub_x_pow_2,-1.0,1.0)
        t_end = time.time()
        print("comp_trep result is ", comp_trep_result, "\tspend time: ", (t_end-t_begin)*1000 / n_test_loops, "ms")

        t_begin = time.time()
        for i in range(0,n_test_loops):
            romberg_result = tester.romberg(normal_functions.exp_x_mul_sqrt_1_sub_x_pow_2,-1.0,1.0)
        t_end = time.time()
        print("romberg rsult is ",romberg_result, "\tspend time: ", (t_end-t_begin)*1000/ n_test_loops, "ms")


        print("#########################\n计算sin(x),在[0,pi/2]区间上的积分\n#######################")

        t_begin = time.time()   
        for i in range(0,n_test_loops):     
            gauss_leg_9_result = tester.gauss_leg_9_common_range(normal_functions.sin_x,0.0,0.5*np.pi)
        t_end = time.time()
        print("gauss_leg_9 result is ", gauss_leg_9_result, "\tspend time: ", (t_end-t_begin)*1000 /n_test_loops, "ms")

        t_begin = time.time()
        for i in range(0,n_test_loops):
            comp_gauss_leg_result = tester.comp_gauss_leg(normal_functions.sin_x,0.0,0.5*np.pi)
        t_end = time.time()
        print("comp_gauss_leg result is ", comp_gauss_leg_result, "\tspend time: ", (t_end-t_begin)*1000 / n_test_loops, "ms")

        t_begin = time.time()
        for i in range(0,n_test_loops):
            comp_trep_result = tester.comp_trep(normal_functions.sin_x,0.0,0.5*np.pi)
        t_end = time.time()
        print("comp_trep result is ", comp_trep_result, "\tspend time: ", (t_end-t_begin)*1000 / n_test_loops, "ms")

        t_begin = time.time()
        for i in range(0,n_test_loops):
            romberg_result = tester.romberg(normal_functions.sin_x,0.0,0.5*np.pi)
        t_end = time.time()
        print("romberg rsult is ",romberg_result, "\tspend time: ", (t_end-t_begin)*1000 / n_test_loops, "ms")


if __name__ == "__main__":
    # s = Integration()
    # s.verbose = True
    # s.gauss_ch1(normal_functions.f_is_one,10)
    # s.gauss_ch2(normal_functions.f_is_one,10)
    # s.gauss_leg_9(normal_functions.f_is_x_power)
    # s.comp_gauss_leg(normal_functions.f_is_x_power,0.,5.)
    #s.comp_trep(normal_functions.sinx_divide_x,0.,1.)
    # result = s.romberg(normal_functions.sinx_divide_x,0.,1.)
    # print("result is ", result)

    tmp_test = Integration()
    tmp_test.verbose = True
    # 细节,关于测试一
    # comp_gauss_leg_result = tmp_test.comp_gauss_leg(normal_functions.exp_x_mul_sqrt_1_sub_x_pow_2,-1.0,1.0)
    # tmp_test.comp_trep(normal_functions.exp_x_mul_sqrt_1_sub_x_pow_2,-1.0,1.0)
    # tmp_test.romberg(normal_functions.exp_x_mul_sqrt_1_sub_x_pow_2,-1.0,1.0)

    # 细节，关于测试二
    # comp_gauss_leg_result = tmp_test.comp_gauss_leg(normal_functions.sin_x,0.0,0.5*np.pi)
    # tmp_test.comp_trep(normal_functions.sin_x,0.0,0.5*np.pi)
    # tmp_test.romberg(normal_functions.sin_x,0.0,0.5*np.pi)


    normal_functions.test_those_methods(1)