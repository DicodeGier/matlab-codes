%%% Test of function

M = load('inst3.mat').M;
K_max = 12;
W_max = 150;
T = 60;
[B,S,k,w] = ARCDP(M,K_max,W_max,T);