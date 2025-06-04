
restart: with(LinearAlgebra): with(plots): Digits := 16:

# Parameters
nbas := 2:  a0 := 6*1: q_vals := [-3.0729830131446736, -0.9270169868553264  ]: p_vals := [0.0, 0.0]:  
# potential
E:= 16; x0 := 1.0; k := 10;
V := E + 0.5*k*(x-x0)^2;
plot_V := plot(V, x = -8.0..10, y = 0..200, color = red, thickness = 2, legend = "Morse Potential");
tay_exp := evalf(convert(taylor(V, x = q_vals, 3), polynom));

# Taylor expansions plots
tay_plots := []:
for q0 in q_vals do
    tay := convert(taylor(V, x = q0, 3), polynom):
    tay_plot := plot(tay, x = -2..4, color = blue, linestyle = dash, thickness = 2, legend = cat("Taylor Expansion around x =", q0)):
    tay_plots := [op(tay_plots), tay_plot]:
end do:
# Original and the Taylor expansion plots
display([plot_V, op(tay_plots)], title = "Original Potential and Taylor Expansions", axes = boxed, thickness = 2);



# Initialize matrices
S := Matrix(nbas, nbas):  
T := Matrix(nbas, nbas):  
Vmat := Matrix(nbas, nbas):  


for i from 1 to nbas do
  for j from i to nbas do
    
    gi := (a0/Pi)^(1/4)*exp(-a0/2*(x - q_vals[i])^2 + I*p_vals[i]*(x - q_vals[i])):
    gj := (a0/Pi)^(1/4)*exp(-a0/2*(x - q_vals[j])^2 + I*p_vals[j]*(x - q_vals[j])):

    
    S[i,j] := evalf(Int(gi*conjugate(gj), x = -infinity .. infinity)):
    S[j,i] := S[i,j]:

    
    T[i,j] := evalf(Int(gi*(-1/2)*diff(conjugate(gj), x$2), x = -infinity .. infinity)):
    T[j,i] := T[i,j]:

    # Potential (LHA around center q_vals)
    V_lha := convert(taylor(V, x = q_vals[i], 3), polynom):
    Vmat[i,j] := evalf(Int(gi * V_lha * conjugate(gj), x = -infinity .. infinity)):
    Vmat[j,i] := Vmat[i,j]:
  end do:
end do:


# Hamiltonian matrix
print(T);
H := T+Vmat;
S:= print(S);
Kinetic_Term := print(T);
Potential_Term := print(Vmat);
Total_Hamiiltonian := print(H);

# Energy:
evals := sort(Eigenvalues(H));
#coefficients
c := Vector[column]([0.010362900006474308+0.0.I,0.010362900006474308+0.0.I]);


E := (Transpose(conjugate(c)).H.c);
#denominator := Re((Transpose(conjugate(c)).S.c));
#E_qtag := numerator / denominator;








