syms x y
log(x) + exp(y)
y(x) = piecewise(x<0, -1, x>0, 1)
clc
syms f(x)
f(x) = x^4-2*x^3+6*x^2-2*x+10
f(-5)
syms y1 y2
y1 = x+3; y2 = 3*x;
solve(y1 == y2)
syms x
solve(x^4 == 1)
syms a b c
subs(cos(a) + sin(b) - exp(2*c), [a b c], [pi/2 pi/4 -1])
syms x f(x)
assume(x>0)
f(x) = 6*x^7-2*x^6+3*x^3-8;
fplot(f)
xlim([-10 10])
ylim([-1e3 1e3])
doubleSol = roots([-8 3 -2 6]) %  double-precision
symsSol = solve(f) % exact. The roots object stores the zeros for symbolic computations
vpaSol = vpasolve(f) % variable-precision
clc
syms x f(x)
assume(x>0)
f(x) = 6*x^7-2*x^6+3*x^3-8;
fplot(f)
xlim([-10 10])
ylim([-1e3 1e3])
doubleSol = roots([-8 3 -2 6]) %  double-precision
symsSol = solve(f) % exact. The roots object stores the zeros for symbolic computations
vpaSol = vpasolve(f) % variable-precision
syms a b y(x)
dsolve(diff(y) == -a*y)
syms a b c d
syms x1 x2
x = [x1; x2];
A = [a b ; c d];
b = A*x
fplot(tan(x))
fplot(cot(x))
fplot(tan(x)+six(x))
fplot(tan(x)+sfix(x))
fplot(tan(x)+sin(x))
fplot(sin(x))
syms t
x = t*sin(5*t);
y = t*cos(5*t);
fplot(x, y)
grid on
syms t
xt = exp(abs(t)/10).*sin(5*abs(t));
yt = exp(abs(t)/10).*cos(5*abs(t));
zt = t;
h = fplot3(xt,yt,zt, [-10,10],'--r');
syms x y
fsurf(sin(x) + cos(y))
fcontour(sin(x) + cos(y))
sym(1/7)
y=sym(y)
y=sym('y')
y=syms(y)
syms y
clear all
A = sym('a', [1 20])
whos
clear all
syms(sym('a', [1 10]))
whos
syms(('a', [1 10]))
whos
syms(('a' [1 10]))
whos
syms f(x,y)
f(x,y) = x^2*y
f(3,2)
dfx = diff(f,x)
dfx(y+1,y)
syms a b c
A = [a b c; c a b; b c a]
sum(A(1,:))
isAlways(sum(A(1,:)) == sum(A(:,2)))
A = sym('A', [2 4])
A = hilb(3)
A = sym(A)
syms x
f = sin(x)^2;
diff(f)
syms x y
f = sin(x)^2 + cos(y)^2;
diff(f)
syms x y
f = sin(x)^2 + cos(y)^2;
diff(f, y)
syms x y
f = sin(x)^2 + cos(y)^2;
diff(f, y, 2)
syms x
f = sin(x)^2;
int(f)
syms x y n
f = x^n + y^n;
int(f)
syms x y n
f = x^n + y^n;
int(f, y)
syms x y n
f = x^n + y^n;
int(f, 1, 10)
syms x
solve(x^3 - 6*x^2 == 6 - 11*x)
syms x y
solve(6*x^2 - 6*x^2*y + x*y^2 - x*y + y^3 - y^2 == 0, y)
syms x y z
[x, y, z] = solve(z == 4*x, x == y, z == x^2 + y^2)
phi = (1 + sqrt(sym(5)))/2;
f = phi^2 - phi - 1
simplify(f)
syms x
f = (x ^2- 1)*(x^4 + x^3 + x^2 + x + 1)*(x^4 - x^3 + x^2 - x + 1);
expand(f)
syms x
g = x^3 + 6*x^2 + 11*x + 6;
factor(g)
syms x
h = x^5 + x^4 + x^3 + x^2 + x;
horner(h)
syms x
f = 2*x^2 - 3*x + 1;
subs(f, 1/3)
syms x
f = x^3 - 15*x^2 - 24*x + 350;
A = [1 2 3; 4 5 6];
subs(f,A)
syms a b c
A = [a b c; c a b; b c a]
alpha = sym('alpha');
beta = sym('beta');
A(2,1) = beta;
A = subs(A,b,alpha)
syms x y
eqn = (x^2 + y^2)^4 == (x^2 - y^2)^2;
fimplicit(eqn, [-1 1])
syms t
fplot3(t^2*sin(10*t), t^2*cos(10*t), t)
syms x y
fsurf(x^2 + y^2)
syms(z)
syms z
assumptions(z)
syms s
f=1/1+s
ilaplace(f)
f=1/(1+s)
ilaplace(f)