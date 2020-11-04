from random import randint

output = ''
BOUND = 20000
N = randint(10, 12500)
output += str(N) + '\n'
for i in range(N):
    first = randint(-BOUND, BOUND)
    second = randint(-BOUND, BOUND)
    third = randint(-BOUND, BOUND)
    fourth = randint(-BOUND, BOUND)
    #print('polygon(({0}, {1}), ({2}, {3}))'.format(first, second, third, fourth))
    output += "{0} {1} {2} {3}\n".format(first, second, third, fourth)

with open("test.txt", "w") as file:
    file.write(output)


