import math
import numpy as np
def equalizeCluster(matrix):
    
    # matrix = [
    #     [394, 378, 397],
    #     [388, 376, 427],
    #     [418, 425, 397]]
    List = np.array(matrix).flatten()
    dim = int(math.sqrt(List.size))
    Nave = int(List.mean())
    misp = sum([abs(i - Nave) for i in List])/2
            
    hist = []
    tot = 0
    for d in range(1,dim):
        it = 1
        curr = []
        while it < List.size:
            index = List.argsort()
            N0, idx0 = List[index[-it]], index[-it]
            if N0 <= Nave:
                it += 1
                continue
            r0,c0 = idx0//dim, idx0%dim
            N1 = float('inf')
            for r in range(dim):
                for c in range(dim):
                    dx = min((r0-r)%dim, (r-r0)%dim)
                    dy = min((c0-c)%dim, (c-c0)%dim)
                    if dx+dy == d and List[r*dim+c] < N1:
                        N1 = List[r*dim+c]
            if N1 >= Nave:
                it += 1
                continue
            idx1 = np.where(List==N1)[0][0]
            dN = min(Nave-N1, N0 - Nave)
            List[idx0] -= dN
            List[idx1] += dN
            tot += dN*d
            curr.append((idx0,idx1,dN))
        hist.append(curr)
    return hist, List, misp,tot

if __name__ == "__main__":
    matrix = [
        [401, 415, 416],
        [393, 385, 391],
        [373, 446, 380]]
    hist, List, misp,tot = equalizeCluster(matrix)
    print(hist)