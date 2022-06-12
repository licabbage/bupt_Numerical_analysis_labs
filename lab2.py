import numpy as np
import matplotlib.pyplot as plt

#######################
# 绘制1/x^2 - c 的图像
######################
def plot(c:int = 2):
        X_axis = np.arange(0.5,2.001,0.001,dtype=np.float32)
        #print(X_axis)
        
        zero = np.zeros(shape= X_axis.shape)
        

        Y_axis = 1/np.power(X_axis,2)-c

        plt.plot(X_axis,Y_axis)
        plt.plot(X_axis,zero)
        plt.show()

        return

if __name__ == "__main__":
    plot()