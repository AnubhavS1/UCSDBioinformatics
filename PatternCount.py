def count(text, sub):
    counter = 0
    for i in range(0, len(text) - len(sub)):
        check = text[i:i + len(sub)]
        if check == sub[0:len(sub)]:
            counter += 1
    return counter


file = open("test.txt", "r")
t = file.readline()
k = file.readline()
print(count(t, k))