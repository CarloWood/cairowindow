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

R = P - Q_final;
error = int(R * R, x);

error_x1_plus_delta = subs(error, x, x1 + delta);
error_x1 = subs(error, x, x1);
def_integral = error_x1_plus_delta - error_x1;
simplified_def_integral = simplify(def_integral);

partial_A = diff(simplified_def_integral, A);
partial_B = diff(simplified_def_integral, B);

eq3 = partial_A == 0;
eq4 = partial_B == 0;

disp(eq3);
