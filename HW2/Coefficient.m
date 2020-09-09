close all;clear;clc
syms a1 a2 b1 b2 g1 g2 z1;
eq1 = g1 + g2 + z1 -1;
eq2 = a1 + a2 + b1 + b2 - 1;
eq3 = b1*(a1 + a2 + b1 + b2) + b2*(a1 + a2 + b2) + a1*a2 -1/2;
eq4 = b2*(g1 + g2 + z1) + g1*(b1 + a2) - 1/2;
eq5 = a1*g2 + b1*g2 -1/2;
eq6 = g1*g2 - 1/2
eq7 = b1 - b2;
eq8 = g1 - g2 + z1;
eqns = [eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8];
solu = solve(eqns, [a1 a2 b1 b2 g1 g2 z1])

a1_s = solu.a1
a2_s = solu.a2
b1_s = solu.b1
b2_s = solu.b2
g1_s = solu.g1
g2_s = solu.g2
z1_s = solu.z1