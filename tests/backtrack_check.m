# Please add `pkg load symbolic;` to your .octaverc file.

# Mark coefficients of used polynomials as global.
global a b c d e f g;

# Declare all symbols that we use.
syms a b c d e f g w0 h;

# Define the parabola P(x) = a + bx + cx^2.
function result = P (x)
  global a b c;
  result = a + b * x + c * x^2;
endfunction

# Define the derivative of P'(x) = b + 2cx.
function result = dP (x)
  global b c;
  result = b + 2 * c * x;
endfunction

# Let w1 be the x-coordinate of the vertex of the P.
w1 = -b / (2 * c);

Lw1 = P(w1) + h;

# Let Q be the cubic Q(x) = d + ex + fx^2 + gx^3, such that
# it goes through (w0, P(w0)), (w1, P(w1)) and Q'(w0) = P'(w0) = b + 2c w0
# and Q''(w0) = P''(w0) = 2c. Then
g = (P(w0) - Lw1 - (w0 - w1) * dP(w0) + (w0 - w1)^2 * c) / (w0 - w1)^3;
f = c - 3 * g * w0;
e = dP(w0) - 2 * f * w0 - 3 * g * w0^2;
d = P(w0) - (e * w0 + f * w0^2 + g * w0^3);

format_str_d = ['d = ', latex(simplify(d))];
disp(format_str_d)
format_str_e = ['e = ', latex(simplify(e))];
disp(format_str_e)
format_str_f = ['f = ', latex(simplify(f))];
disp(format_str_f)
format_str_g = ['g = ', latex(simplify(g / h))];
disp(format_str_g)


# Calculate the values of g, f, e and d when Lw1 equals the y-coordinate of the vertex of P.
g_match = subs(g, Lw1, P(w1));
g_match = simplify(g_match);
f_match = subs(f, Lw1, P(w1));
f_match = simplify(f_match);
e_match = subs(e, Lw1, P(w1));
e_match = simplify(e_match);
d_match = subs(d, Lw1, P(w1));
d_match = simplify(d_match);

# Print the result.
g_match_str = char(g_match);
f_match_str = char(f_match);
e_match_str = char(e_match);
d_match_str = char(d_match);
format_str = ['When L(w1) equals the y-coordinate of the vertex, then g = ', g_match_str, ', f = ', f_match_str, ', e = ', e_match_str, ' and d = ', d_match_str, '.'];
disp(format_str)

# Define Lw1 relative to P(w1).
%x = sym('x');

%zero = (-f + sqrt(f^2 - 3*g*e)) / (3*g);
%zero_x = subs(zero, Lw1, P(w1) + x);
%zero_x = simplify(zero_x);

%disp(zero_x)

%zero_x2 = simplify(expand(zero_x))
%disp(zero_x2)

