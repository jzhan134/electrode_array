import random, equalizeCluster, numpy as np, matplotlib.pyplot as plt
d = 5
misp_hist = []
move_hist = []
for it in range(10000):
    if d == 3:
        candidate = list(range(179*179))
        sample = random.sample(candidate,3600)
        cnt = [[0 for _ in range(d)] for _ in range(d)]
        for it in sample:
            r,c = it//180 - 90, it%180 - 90
            cnt[(r+30)//60+1][(c+30)//60+1] += 1
    elif d == 4:
        candidate = list(range(159*159))
        sample = random.sample(candidate,3600)
        cnt = [[0 for _ in range(d)] for _ in range(d)]
        for it in sample:
            r,c = it//160 - 80, it%160 - 80
            cnt[r//40+2][c//40+2] += 1
    else:
        candidate = list(range(199*199))
        sample = random.sample(candidate,3600)
        cnt = [[0 for _ in range(d)] for _ in range(d)]
        for it in sample:
            r,c = it//200-100, it%200-100
            cnt[(r+20)//40+1][(c+20)//40+1] += 1

    hist, List,misp, tot = equalizeCluster.equalizeCluster(cnt)
    misp_hist.append(misp)
    move_hist.append(tot)
    # print(np.array(cnt).reshape((d,d)))
    # print(hist)
    # print(misp,tot)
file1 = open("d5.txt","w") 
file1.writelines(' '.join([str(i) for i in misp_hist])+'\n')
file1.writelines(' '.join([str(i) for i in move_hist]))
file1.close() 