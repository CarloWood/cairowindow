syms x x1 a b c d e A B C D delta;

disp("Let P(x) be a fourth degree polynomial:");
P = a * x^4 + b * x^3 + c * x^2 + d * x + e;
disp(P);

disp("\nLet Q(x) be a third degree polynomial:");
Q = A * x^3 + B * x^2 + C * x + D;
disp(Q);

disp("\nThen P'(x) = ");
dP = diff(P, x);
disp(dP);
disp("\nand Q'(x) = ");
dQ = diff(Q, x);
disp(dQ);

disp("\nBy setting P(x1) = Q(x1), and P'(x1) = Q'(x1):");
eq1 = subs(P, x, x1) == subs(Q, x, x1);
eq2 = subs(dP, x, x1) == subs(dQ, x, x1);
disp(eq1);
disp(eq2);

disp("\nwe can solve for C and D to get:");
C_sol = solve(eq2, C);
eq1_subC = subs(eq1, C, C_sol);
D_sol = solve(eq1_subC, D);
Q_final = subs(Q, {C, D}, {C_sol, D_sol});
disp(Q_final);

disp("\n1260 times the integral from x1 to x1+δ of (P(x) - Q(x))² is given by:");
R = P - Q_final;
error = int(R * R, x);
error_x1_plus_delta = subs(error, x, x1 + delta);
error_x1 = subs(error, x, x1);
def_integral = error_x1_plus_delta - error_x1;
simplified_def_integral = simplify(def_integral);
disp(simplify(1260 * simplified_def_integral));

partial_A = diff(simplified_def_integral, A);
partial_B = diff(simplified_def_integral, B);

disp("\nSetting the partial derivative with respect to A to zero gives:");
eq3 = partial_A == 0;
disp(eq3);

disp("\nSetting the partial derivative with respect to B to zero gives:");
eq4 = partial_B == 0;
disp(eq4);

disp("\nFrom this we can solve A and B:");
B_sol = solve(eq4, B);
eq3_subB = subs(eq3, B, B_sol);
A_sol = solve(eq3_subB, A);
B_sol_subA = subs(B_sol, A, A_sol);
B_sol_subA_simplified = simplify(B_sol_subA);
disp("\nA =");
disp(A_sol);
disp("\nB =");
disp(B_sol_subA_simplified)

ddQ = diff(dQ, x);
ddQ_x1 = subs(ddQ, x, x1);
ddQ_x1_subAB = simplify(subs(subs(ddQ_x1, A, A_sol), B, B_sol_subA));
dddQ = diff(ddQ, x);
dddQ_x1_subAB = simplify(subs(dddQ, A, A_sol));
ddP = diff(dP, x);
dddP = diff(ddP, x);
disp("\nQ''(x1) = ");
disp(ddQ_x1_subAB)
disp("\nQ'''(x1) = ");
disp(dddQ_x1_subAB)

disp("\nNote that P''(x1) =");
ddP_x1 = simplify(subs(ddP, x, x1));
disp(ddP_x1);

disp("\nAnd P'''(x1) =");
dddP_x1 = simplify(subs(dddP, x, x1));
disp(dddP_x1);
