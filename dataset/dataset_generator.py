import random
def generator(length1,length2,filename1,filename2):
    pool = ["A","T","C","G"]
    seq1 = []
    seq2 = []
    for i in range(length1):
        seq1.append(random.choice(pool))
    for i in range(length2):
        seq2.append(random.choice(pool))
    file = open(filename1,"w")
    for c in seq1:
        file.write(c)
    file = open(filename2, "w")
    for c in seq2:
        file.write(c)
if __name__ == "__main__":
    generator(20000,20000,"seq1.txt","seq2.txt")
