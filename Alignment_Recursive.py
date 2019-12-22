import numpy as np
import os
import time
import sys
import psutil

#delta fitting is the scoring function

keys = ['A', 'C', 'T', 'G']
delta_fitting = {}
for i in range(len(keys)):
    delta_fitting[keys[i]] = {k : v for (k,v) in zip(keys, [3 if keys[i] == keys[j]  else -1 for j in range(len(keys))])}
delta_fitting['-'] = {}
delta_fitting['-']['-'] = -2
for key in keys:
    delta_fitting[key]["-"] = -2
    delta_fitting["-"][key] = -2

def score(x,y):
    prev_row = [(i * delta_fitting["-"]["-"]) for i in range(len(y) + 1)]
    cur_row = [float("-inf")] * (len(y) + 1)
    for i in range(1, len(x) + 1):
        for j in range(0, len(y) + 1):
            if j > 0:
                cur_row[j] = prev_row[j - 1] + delta_fitting[x[i - 1]][y[j - 1]]

                cur_row[j] = max(cur_row[j], cur_row[j - 1] + delta_fitting["-"][y[j - 1]])
            cur_row[j] = max(cur_row[j], prev_row[j] + delta_fitting[x[i - 1]]["-"])

        prev_row = cur_row
        cur_row = [float("-inf")] * (len(y) + 1)
    return prev_row

def Hirschberg(x,y):
    z = []
    w = []
    if not x:
        for i in range(len(y)):
            z.append("-")
            w.append(y[i])
    elif not y:
        for i in range(len(x)):
            z.append(x[i])
            w.append("-")
    elif len(x) == 1 or len(y) == 1:
        return NeedlemanWunsch(x,y,delta_fitting)
    else:
        xlen = len(x)
        xmid = len(x) // 2
        ylen = len(y)

        ScoreL = np.array(score(x[:xmid], y))
        ScoreR = score(x[xmid:xlen][::-1], y[::-1])
        ScoreR.reverse()
        ScoreR = np.array(ScoreR)
        ymid = np.argmax(ScoreR + ScoreL)
        z1,w1 = Hirschberg(x[:xmid], y[:ymid])
        z2,w2 = Hirschberg(x[xmid: xlen], y[ymid: ylen])
        z += z1
        z += z2
        w += w1
        w += w2
        return z,w
    return z,w


UP = (-1,0)
LEFT = (0, -1)
TOPLEFT = (-1, -1)
ORIGIN = (0, 0)


def traceback_global(v, w, pointers):
    i, j = len(v), len(w)
    new_v = []
    new_w = []
    while True:
        di, dj = pointers[i][j]
        if (di, dj) == LEFT:
            new_v.append('-')
            new_w.append(w[j - 1])
        elif (di, dj) == UP:
            new_v.append(v[i - 1])
            new_w.append('-')
        elif (di, dj) == TOPLEFT:
            new_v.append(v[i - 1])
            new_w.append(w[j - 1])
        i, j = i + di, j + dj
        if i <= 0 and j <= 0:
            break

    return new_v[::-1], new_w[::-1]


def NeedlemanWunsch(v, w, delta):

    if len(v) < len(w):
        v,w = w,v
    M = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    pointers = [[ORIGIN for j in range(len(w) + 1)] for i in range(len(v) + 1)]

    for i in range(1, len(w) + 1):
        M[0][i] += M[0][i - 1] + delta[w[i - 1]]['-']
        pointers[0][i] = LEFT
    for i in range(1, len(v) + 1):
        M[i][0] += M[i - 1][0] + delta[v[i - 1]]['-']
        pointers[i][0] = UP

    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            M[i][j] = -999999
            if M[i - 1][j] + delta[v[i - 1]]['-'] > M[i][j]:
                M[i][j] = M[i - 1][j] + delta[v[i - 1]]['-']
                pointers[i][j] = UP
            if M[i][j - 1] + delta[w[j - 1]]['-'] > M[i][j]:
                M[i][j] = M[i][j - 1] + delta[w[j - 1]]['-']
                pointers[i][j] = LEFT
            if delta[v[i - 1]][w[j - 1]] + M[i - 1][j - 1] > M[i][j]:
                M[i][j] = delta[v[i - 1]][w[j - 1]] + M[i - 1][j - 1]
                pointers[i][j] = TOPLEFT


    alignment = traceback_global(v, w, pointers)
    return alignment

def Hirschberg_Entrance(v,w):
    if len(v) < len(w):
        v,w = w,v
    return Hirschberg(v,w)
def testcase(filename1="dataset/seq1.txt",filename2="dataset/seq2.txt"):
    seq1 = open(filename1,"r").readline()
    seq2 = open(filename2,"r").readline()
    return seq1,seq2

def evaluation(st1,st2,delta):
    total = 0

    for c1,c2 in zip(st1,st2):
        total += delta[c1][c2]
    return total

if __name__ == "__main__":
    process = psutil.Process(os.getpid())

    filename1, filename2 = sys.argv[1], sys.argv[2]
    v, w = testcase(filename1="dataset/seq1_size{}.txt".format(filename1),
                    filename2="dataset/seq2_size{}.txt".format(filename2))

    print("Following is the performance test for size {} * {}".format(filename1, filename2))
    st1 = time.time()
    ali1, ali2 = Hirschberg_Entrance(v, w)
    st2 = time.time()
    print("Hirschberg Performance:")
    print("time usage: {} s".format(st2 - st1))
    print("memory usage: {} MB(unit might change with system, MB for mac OS)".format(process.memory_info().rss / (1024 * 1024)))

    st1 = time.time()
    cali1, cali2 = NeedlemanWunsch(v, w, delta_fitting)
    st2 = time.time()
    print("NeedlemanWunsch Performance:")
    print("time usage: {} s".format(st2 - st1))

    print("memory usage: {} MB(unit might change with system, MB for mac OS)".format(
        process.memory_info().rss / (1024 * 1024)))
    if evaluation("".join(ali1), "".join(ali2), delta_fitting) == evaluation("".join(cali1), "".join(cali2),delta_fitting) and len(ali1) == len(cali1):

        print("Hirschberg returned optimal alignment")
    else:
        print("Hirschberg returned an alignment which has a score different from the one by NeedlemanWunsch")
    # print("".join(ali1))
    print("Finished the performance test for size {} * {}".format(filename1, filename2))
    print("\n\n")
