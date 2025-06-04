
restart; with(plots): with(LinearAlgebra): Digits := 16:
psi := (2*alpha/Pi)^(1/4)*exp(-alpha*(x)^2):
rho:=simplify(psi*psi);
t10:=solve(rho=rho_cutoff, x);
alpha := 9.16:
S := 0.95:

# Parameters
a := [9*alpha]:          
S := 0.92:                 
n := nops(a):
rho_cutoff := 1e-1;

delta_x := t10[1];
x_min := -delta_x:
x_max := delta_x:

for i from 1 to n do
    al := a[i]: 
    dx := sqrt(-4/al*ln(S)):
    print(dx);
    N := round((x_max - x_min)/dx)+1:
end do;

# Density;
x0:=0; alpha := 9.16;
psi := (alpha/Pi)^(1/4)*exp(-alpha*(x-x0)^2/2):
rho:=simplify(psi*psi);
solve(rho=rho_cutoff, x);
delta_x := sqrt((1/alpha)*ln(sqrt(alpha/Pi)/rho_cutoff));
cutoff_line := rho_cutoff;



# Density
p1 := plot(rho, x = -1.5..1.5, color = blue, thickness = 2, legend = "rho(x)"):
p2 := plot(cutoff_line, x = -1.5..1.5, color = red, linestyle = dash, legend = "cutoff"):
p3 := plot([[delta_x, 0], [delta_x, eval(rho, x = delta_x)]], color = green, linestyle = dash, thickness = 2):
p4 := plot([[-delta_x, 0], [-delta_x, eval(rho, x = -delta_x)]], color = green, linestyle = dash, thickness = 2):
print(delta_x);
print(-delta_x);

display([p1, p2, p3, p4], labels = ["x", "rho(x)"]);



