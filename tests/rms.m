syms x Qx1 dQx1 ddQx1 dddQx1 u x1;

disp("\nLet Q(x) be a third degree polynomial:");
Q = Qx1 + dQx1 * (x - x1) + (1/2) * ddQx1 * (x - x1)^2 + (1/6) * dddQx1 * (x - x1)^3;
disp(Q);

disp("\nQ'(x) = ");
dQ = diff(Q, x);
disp(simplify(dQ));

disp("\nthe integral from x1 to x2 of Q'(x)² divided by (x2 - x1) is given by (where u = x2 - x1):");
integral_square = int(dQ^2, x);
integral_squared_x1 = subs(integral_square, x, x1);
integral_squared_x2 = subs(integral_square, x, u + x1);
avg = (integral_squared_x2 - integral_squared_x1) / u;
simplified_avg = simplify(avg);
disp(simplified_avg);
