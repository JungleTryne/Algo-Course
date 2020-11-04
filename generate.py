import os
from random import randint

for j in range(1000):
    output = ''
    BOUND = 1000
    N = randint(1, 5)
    output += str(N) + '\n'
    for i in range(N):
        first = randint(-BOUND, BOUND)
        second = randint(-BOUND, BOUND)
        third = randint(-BOUND, BOUND)
        fourth = randint(-BOUND, BOUND)
        output += "{0} {1} {2} {3}\n".format(first, second, third, fourth)

    with open("test.txt", "w") as file:
        file.write(output)
    os.system("cmake-build-debug/module3 < test.txt")
    print(j)