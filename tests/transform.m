syms c0 c1 c2 c3 x xs xt;

% Test with (x - 37) * (x + 13) * (x - 127):
%
%           3        2                 
%   P(x) = x  - 151⋅x  + 2567⋅x + 61087
c0_val = vpa(61087);
c1_val = vpa(2567);
c2_val = vpa(-151);
c3_val = vpa(1);

disp("Let p(x) (or px) be the cubic:");
px = c0 + c1 * x + c2 * x^2 + c3 * x^3;
disp(px);
px_substituted = subs(px, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val]);
disp(["\nTest case: P(x) = ", char(vpa(px_substituted))]);

disp("\nSubstitute x = xs - c2 / (3 c3) to get f(xs) =");
fxs = collect(expand(simplify(subs(px, x, xs - c2 / (3 * c3)))), xs);
disp(fxs);
x_subs = xs - c2 / (3 * c3);
Ix = vpa(subs(-c2 / (3 * c3), [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val]));
disp(["\nSubstitute: x --> ", char(vpa(subs(x_subs, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val])))]);
fxs_substituted = subs(fxs, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val]);
disp(["\nTest case: f(xs) = ", char(vpa(fxs_substituted))]);

disp("\nThis function has an extreme k at");
dfxs = diff(fxs, xs);
extrema_eq = dfxs == 0;
extrema = solve(extrema_eq, xs);
k = extrema(2);
disp(k);
k_substituted = subs(k, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val]);
disp(["\nTest case: k = ", char(vpa(k_substituted))]);

disp("\nSubstitute xs = k * xt to get g(xt) =");
gxt = collect(expand(simplify(subs(fxs, xs, xt * k))), xt);
disp(gxt);
xs_subs = k * xt;
disp(["\nSubstitute: xs --> ", char(vpa(subs(xs_subs, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val])))]);
gxt_substituted = subs(gxt, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val]);
disp(["\nTest case: g(xt) = ", char(vpa(gxt_substituted))]);

disp("\nThe coefficient of cubic term (gxtc3) is:");
dgxt = diff(gxt, xt);
ddgxt = diff(dgxt, xt);
dddgxt = diff(ddgxt, xt);
gxtc3 = dddgxt / 6;
disp(gxtc3);

disp("\nDivide g(xt) by this coefficient; h(xt) = g(xt) / gxtc3, which should give a monic polynomial:");
hxt = collect(expand(simplify(gxt / gxtc3)), xt);
dhxt = diff(hxt, xt);
ddhxt = diff(dhxt, xt);
dddhxt = diff(ddhxt, xt);
hxtc3 = dddhxt / 6;
disp("coefficient of tx^3 (should be one):");
disp(simplify(hxtc3));

hxt2 = hxt - hxtc3 * xt^3;
dhxt2 = diff(hxt2, xt);
ddhxt2 = diff(dhxt2, xt);
disp("coefficient of tx^2 (should be zero):");
disp(ddhxt2);
disp("coefficient of tx (should be -3):");
disp(simplify(dhxt2));

disp("\nThe constant term of h(xt) is:");
hxtc0 = simplify(subs(hxt, xt, 0));
disp(hxtc0);

hxt_substituted = subs(hxt, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val]);
disp(["\nTest case: h(xt) = ", char(vpa(hxt_substituted))]);

% Test with (x - 37) * (x + 13) * (x - 127):
%
%           3        2                 
%   P(x) = x  - 151⋅x  + 2567⋅x + 61087
c0_val = vpa(61087);
c1_val = vpa(2567);
c2_val = vpa(-151);

digits(32);
%c0_val = vpa(103503) / vpa(10000);
%c1_val = vpa(99721) / vpa(1000);
%c2_val = vpa(172963) / vpa(10000);
hxtc0_substituted = subs(hxtc0, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val]);
result = vpa(hxtc0_substituted);
disp(["\nThe result of the expression is: ", char(result)]);

disp("\nThe transformed cubic is therefore:");
px_transformed = result - 3 * x + x^3;
disp(px_transformed);

disp("\nwhich has roots:");
eq = px_transformed == 0;
roots = real(solve(eq, x));
disp(roots);

disp(["\nMultiply with ", char(k_substituted)]);
roots_mult = k_substituted * roots;
disp(roots_mult);

disp(["\nAdd back ", char(Ix)]);
roots_done = roots_mult + Ix;
for i = 1:length(roots_done)
  printf("%.6f\n", double(roots_done(i)));
end

disp("\nVictor!");

% Next, test with (-3652139/3) + (37901/3)⋅x - 151⋅x² + x³
%
c0_val = sym(-3652139)/3;
c1_val = sym(37901)/3;
c2_val = sym(-151);
c3_val = sym(1);

digits(8);

px_substituted = subs(px, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val]);
disp("\nNext test case: P(x) =\n");
disp(px_substituted);

Ix = subs(-c2 / (3 * c3), [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val]);
disp(["\nSubstitute: x --> ", char(vpa(subs(x_subs, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val])))]);
fxs_substituted = subs(fxs, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val]);
disp(["\nTest case: f(xs) = ", char(vpa(fxs_substituted))]);

% Extract C0 and C1 from fxs = C0 + C1 xs + c3 xs^3
% dfxs = C1 + 3 c3 xs^2
disp("\nLet k =");
C0 = simplify(subs(fxs, xs, 0));
C1 = dfxs - 3 * c3 * xs^2;
k = sqrt(C1 / 3);
disp(k);
k_substituted = subs(k, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val]);
disp(["\nTest case: k = ", char(vpa(k_substituted))]);
D0 = C0 / k^3;
disp("\nLet h(xt) =");
hxt = D0 + 3 * xt + xt^3
hxt_substituted = subs(hxt, [c0, c1, c2, c3], [c0_val, c1_val, c2_val, c3_val]);
disp(["\nTest case: h(xt) = ", char(vpa(hxt_substituted))]);

D0_div_c3 = subs(D0, c3, 1);
disp("\nFor the monic case the constant term is:");
disp(D0_div_c3);
