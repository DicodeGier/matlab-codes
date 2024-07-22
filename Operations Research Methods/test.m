clear
clc

parent1 = [1 3 5 2 4 6]
parent2 = [6 3 2 5 4 1]

all_numbers = 1:6

offspring1 = zeros(1,6)
offspring1(1:3) = parent1(1:3)
to_find = all_numbers(~ismember(1:6,parent1(1:3)))
offspring1(4:6) = intersect(parent2,to_find,'stable')

offspring2 = zeros(1,6)
offspring2(1:3) = parent2(1:3)
to_find2 = all_numbers(~ismember(1:6,parent2(1:3)))
offspring2(4:6) = intersect(parent1,to_find2,'stable')
