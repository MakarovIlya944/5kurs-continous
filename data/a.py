def qq():
    a = [[0,0,0],[1,0,0],[0,1,0],[1,1,0],
    [1,1,1],[1,1,0],[1,1,0],[1,1,0]]

    b = [[1,1,1],
    [1,1,2],
    [1,2,1],
    [2,1,1],
    [1,2,2],
    [2,1,2],
    [2,2,1],
    [2,2,2]]

    c = [
    [0,0,0],
    [2,0,0],
    [0,2,0],
    [2,2,0],
    [0,0,2],
    [2,0,2],
    [0,2,2],
    [2,2,2]]

    with open('out','w') as f:
        for i in c:
            f.writelines([f'{i[0]+j[0]}{i[1]+j[1]}{i[2]+j[2]}\n' for j in b])
            # f.writelines(['\n'])

def u(i):
    t = i // 4
    t = t % 2
    t *= 2
    t += i % 2
    return t + 1

def v(i):
    t = i // 8
    t *= 2
    t += (i // 2) % 2
    return t + 1

for i in range(16):
    t = u(i)
    t = v(i)
    print(str(u(i)) + str(v(i)))