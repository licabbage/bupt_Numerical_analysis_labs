import numpy as np
import struct
def f(input:float) -> float:
    return input*2


def callf(f,input:float):
    print (f(input))
    return

if __name__ == "__main__":
    x = np.ones(5,dtype=np.float32)
    y = np.array([2,3,4,5,6],dtype=np.float32)
    print(y)


    x[0] = 1
    x[2] = 3

    # print(np.power(x,2))

    print(x)
    
    print(x*y)
    np.multiply(x,y).sum

    curr_n = int(10)
    curr_n *=2
    assert(curr_n == x.size * 4)
    print(x.size)


    print(16**0.5)

    print(0.9460833108884718 - 0.9456908635827013)

    s = float(0.125)
    s = np.float32(0.125)
    print(struct.pack(">f",s))

