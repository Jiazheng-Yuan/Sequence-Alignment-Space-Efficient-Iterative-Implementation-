import numpy as np

# keys = ['A', 'C', 'T', 'G', '-']
# delta_fitting = {}
# for i in range(len(keys)):
#     delta_fitting[keys[i]] = {k : v for (k,v) in zip(keys, [1 if keys[i] == keys[j]  else -1 for j in range(len(keys))])}
# delta_fitting['-']['-'] = -1
keys = ['A', 'C', 'T', 'G']
delta_fitting = {}
for i in range(len(keys)):
    delta_fitting[keys[i]] = {k : v for (k,v) in zip(keys, [3 if keys[i] == keys[j]  else -1 for j in range(len(keys))])}
delta_fitting['-'] = {}
delta_fitting['-']['-'] = -2
for key in keys:
    delta_fitting[key]["-"] = -2
    delta_fitting["-"][key] = -2
def score(x,y,match=1,mismatch=-1):
    # if len(x) > len(y):
    #     x,y = y,x
    # prev_col = [(i * mismatch) for i in range(len(x) + 1)]
    # cur_col = [float("-inf")] * (len(x) + 1)
    #
    # for j in range(1,len(y) + 1):
    #     for i in range(0,len(x) + 1):
    #         if i > 0:
    #             cur_col[i] = prev_col[i - 1]
    #             cur_col[i] += match if x[i - 1] == y[j - 1] else mismatch
    #             cur_col[i] = max(cur_col[i], cur_col[i - 1] + mismatch)
    #         cur_col[i] = max(cur_col[i],prev_col[i] + mismatch)
    #
    #     prev_col = cur_col
    #     cur_col = [0] * (len(x) + 1)
    #
    # return prev_col

    # prev_row = [(i * mismatch) for i in range(len(y) + 1)]
    # cur_row = [float("-inf")] * (len(y) + 1)
    # for i in range(1,len(x) + 1):
    #     for j in range(0,len(y) + 1):
    #         if j > 0:
    #             cur_row[j] = prev_row[j - 1]
    #             cur_row[j] += match if x[i - 1] == y[j - 1] else -1
    #             cur_row[j] = max(cur_row[j], cur_row[j - 1] + mismatch)
    #         cur_row[j] = max(cur_row[j],prev_row[j] + mismatch)
    #
    #     prev_row = cur_row
    #     row_row = [float("-inf")] * (len(y) + 1)
    #
    # return prev_row
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
    id = 0
    if len(x) < len(y):
        x,y = y,x
    result = {}
    visited = {}
    #row_id,col_id = 0,0

    stack = [(x,y,id)]
    original = (x,y,id)
    while stack:
        top = stack.pop()
        x,y,useless = top
        z = []
        w = []
        if top not in visited:
            id += 1
            if not x:
                for i in range(len(y)):
                    z.append("-")
                    w.append(y[i])
                result[top] = [z,w]

            elif not y:
                for i in range(len(x)):
                    z.append(x[i])
                    w.append("-")
                result[top] = [z, w]
            elif len(x) == 1 or len(y) == 1:
                z,w = NeedlemanWunsch(x,y,delta_fitting)
                result[top] = [z, w]
            else:
                xlen = len(x)
                xmid = len(x) // 2
                ylen = len(y)

                ScoreL = np.array(score(x[:xmid], y))
                ScoreR = score(x[xmid:xlen][::-1], y[::-1])
                ScoreR.reverse()
                ScoreR = np.array(ScoreR)
                ymid = np.argmax(ScoreR + ScoreL)

                stack.append(top)
                x1,x2 = x[:xmid],x[xmid: xlen]
                y1,y2 = y[:ymid],y[ymid: ylen]
                stack.append((x1, y1, id))
                id += 1
                stack.append((x2, y2, id))
                id += 1
                visited[top] = [(x1, y1,id - 2),(x2, y2,id - 1)]
        else:
            z1, w1 = result[visited[top][0]]
            z2, w2 = result[visited[top][1]]
            z += z1
            z += z2
            w += w1
            w += w2
            result[top] = [z,w]
            result.pop(visited[top][0])
            #if visited[top][1] in result:
            result.pop(visited[top][1])
            visited.pop(top)

    return result[original]


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
        if (i <= 0 and j <= 0):
            break

    return new_v[::-1], new_w[::-1]


def NeedlemanWunsch(v, w, delta):
    """
    Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment as
    computed by traceback_global.

    :param: v
    :param: w
    :param: delta
    """
    M = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    pointers = [[ORIGIN for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    score, alignment = None, None
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

    # YOUR CODE HERE

    alignment = traceback_global(v, w, pointers)
    score = M[-1][-1]
    #     for row in M:
    #         print(M)
    return alignment
def testcase(filename1="dataset/seq1.txt",filename2="dataset/seq2.txt"):
    seq1 = open(filename1,"r").readline()
    seq2 = open(filename2,"r").readline()
    return seq1,seq2
import os

def evaluation(st1,st2,delta):
    total = 0
    for c1,c2 in zip(st1,st2):
        total += delta[c1][c2]
    return total
import time,resource
if __name__ == "__main__":
    #os.system("ls")
    #v,w = "ATCGAAAAATCG","TATGGGGGATG"
    #v,w = testcase(filename1="Sequence-Alignment-Space-Efficient-Iterative-Implementation-/dataset/seq1.txt",filename2="Sequence-Alignment-Space-Efficient-Iterative-Implementation-/dataset/seq2.txt")
    v,w = testcase()
    st1 = time.time()
    ali1,ali2 = Hirschberg(v,w)
    st2 = time.time()
    print(st2 - st1)
    print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024))
    st1 = time.time()
    cali1,cali2 = NeedlemanWunsch(v,w,delta_fitting)
    st2 = time.time()
    print(st2 - st1)
    print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024))
    print(evaluation("".join(ali1),"".join(ali2),delta_fitting) == evaluation("".join(cali1),"".join(cali2),delta_fitting))
    # print("".join(ali1))
    # print("".join(ali2))
    # print("".join(cali1))
    # print("".join(cali2))
    #print(score("aaaaa","baa"))