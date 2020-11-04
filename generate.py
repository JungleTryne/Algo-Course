from random import randint

output = ''
BOUND = 1000
N = randint(1, 5)
output += str(N) + '\n'
for i in range(N):
    first = randint(-BOUND, BOUND)
    second = randint(-BOUND, BOUND)
    third = randint(-BOUND, BOUND)
    fourth = randint(-BOUND, BOUND)
    #print('polygon(({0}, {1}), ({2}, {3}))'.format(first, second, third, fourth))
    output += "{0} {1} {2} {3}\n".format(first, second, third, fourth)

#output += "50 50 51 51"

with open("test.txt", "w") as file:
    file.write(output)


