# bupt_Numerical_analysis_labs
北邮数值计算与符号分析实验

北邮数值计算与符号分析期中加期末实验
 
lab1.py 是期中实验问题甲：
令
    x0 = (pc + q) / (c + r),
求  
    ((x0-sqrt(c)) /(x0 + sqrt(c)) )^2
的最小无穷范数,给出p,q,r的值以及最小无穷范数结果

record.txt 是期中实验的一些中间输出结果
cont_2_bits.py 用来数二进制0的个数，无特殊意义，可忽略

lab2.py是期末实验问题乙，仅仅绘制了一个图片用于分析

integration_methods.py 是期末实验：
    包括了要求实现的几种积分方法，包括
        1）高斯-切比雪夫 Ⅰ型
        2）高斯-切比雪夫 Ⅱ型
        3.1）高斯-勒让德（gause-legendre) 9点,[-1,1]区间
        3.2)高斯-勒让德（gause-legendre) 9点,任意区间
        4.1)一般区间的2点 高斯-勒让德（gause-legendre） 公式  func为函数,[a,b]为一般区间
        4.2)复化的 逐次减半 2点高斯 高斯-勒让德（gause-legendre） 公式
        5)逐次减半复化梯形公式
        6)Romberg求积公式

以及用这些方法去测试积分结果
    六种方式测试 f(x) = exp(x) * sqrt(1-x^2) 在区间[-1,1]的积分结果
    3，4，5，6种方式测试f(x) =sin(x) 在区间[0,pi/2]的积分结果

test.py为因使用python不熟练所测试一些函数用的代码，可忽略


