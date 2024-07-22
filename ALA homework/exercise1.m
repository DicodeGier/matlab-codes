clear
n = input('number of values ')
vector_numbers = []
for i=1:n
    userResponse = input('enter the number ')
    vector_numbers(i) = userResponse
end

equation_vector = [syms 'y_k%d' 'y_k%d' [1 2]]


