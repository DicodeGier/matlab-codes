clear
clc

addpath("C:\gurobi1001\win64\matlab")

model = []
model.obj = [1 1]
model.A = sparse([1,1])
model.rhs = max([1,2])
model.sense = '<'
model.modelsense = 'max'
model.vtype = ['C','C']
solution = gurobi(model)
