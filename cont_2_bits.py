
if __name__ == "__main__":
    tmp = 0.0000000000000000000000001000
    counter = int(0)
    while tmp < 1.0:
        tmp *=10
        counter+=1

    print(counter)