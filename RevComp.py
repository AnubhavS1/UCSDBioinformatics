def reverseComplement (strand):
    complement = ""
    for i in range(len(strand)):
        if strand[i] == "A":
            complement += "T"
        elif strand[i] == "T":
            complement += "A"
        elif strand[i] == "C":
            complement += "G"
        elif strand[i] == "G":
            complement += "C"
    recComp = ""
    for i in range(len(complement)):
        recComp += complement[len(complement) - i - 1]
    return recComp

file = open("dataset_3_2.txt", "r")
s = file.readline()
print(s)
print(reverseComplement(s))
