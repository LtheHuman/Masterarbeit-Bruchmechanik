from autofreefem import *
from autofreefem.utils import traverse_freefem, rep_derivatives
from sympy import latex

name = 'plate'
x, y = sp.symbols('x y')
f = Constant(r'f', 18.27 * 1e-3, 'load', kilogram/(meter*second))
EI = Constant(r'EI', 1.0, 'EI', kilogram*meter**2)
D = Constant(r'D', 0.1, 'D', meter)
Th = Domain(r'\Omega','Th')

# fluid pressure
u = ScalarField(r'u', P2, Th, '0.', [])
v = ScalarField(r'\delta u', P2, Th)
#
gradu = inner(identity()-TensorProduct(SurfaceNormalVector(),SurfaceNormalVector()),grad(u))
gradgradu = inner(identity()-TensorProduct(SurfaceNormalVector(),SurfaceNormalVector()),grad(gradu))
gradv = inner(identity()-TensorProduct(SurfaceNormalVector(),SurfaceNormalVector()),grad(v))
gradgradv = inner(identity()-TensorProduct(SurfaceNormalVector(),SurfaceNormalVector()),grad(gradv))

# only plate
lag = Lagrangian([u], [v], dx(f*u, Th),
      dx(D*inner2(grad(grad(v)), grad(grad(u))), Th) - dx(f*v, Th)
     , dimensions=2)
lag.compute()

#lag.setup_optimization('plate', 'mesh.txt',boundaryLabels=[1],fixedLabels=[2])

exp = sp.simplify(
    (lag.shape_Derivative
                  .subs(lag.varAdjoint[0],-lag.varAnsatz[0])
                  .subs(SurfaceNormalVector(),ex()) ) )

print(latex(exp))

s = traverse_freefem(exp).simplify()
print("plate only:",rep_derivatives(str(s.subs(sp.Function('Vy')(x, y),0).subs(sp.Function('Vx')(x, y),sp.Function('Vx')(x)).doit())
                                    ,lag.fields).replace(f'Derivative(Vx({x}), x)', 'dx(Vx)').replace(f'Derivative(Vx({x}), (x, 2))', 'dxx(Vx)')) # for vertical stiffener          

# only stiffener vertical
lag = Lagrangian([u], [v], 0, dsx(EI*inner2(gradgradu,gradgradv), Th, 2) ,                   
                 dimensions=2)
lag.compute()
exp = sp.simplify(
    (lag.shape_Derivative
                  .subs(lag.varAdjoint[0],-lag.varAnsatz[0])
                  .subs(SurfaceNormalVector(),ex()) ) )
from autofreefem.utils import traverse_freefem
x, y = sp.symbols('x y')
s = traverse_freefem(exp).simplify()
print("only stiffener vertical:",s.subs(sp.Function('Vy')(x, y),0).subs(sp.Function('Vx')(x, y),sp.Function('Vx')(x)).doit())

# only sensitivity for a horizontal stiffener for movement of a vertical stiffener
lag = Lagrangian([u], [v], 0, dsx(EI*inner2(gradgradu,gradgradv), Th, 2),                   
                 dimensions=2)
lag.compute()
exp = sp.simplify(
    (lag.shape_Derivative
                  .subs(lag.varAdjoint[0],-lag.varAnsatz[0])
                  .subs(SurfaceNormalVector(),ey()) ) )

x, y = sp.symbols('x y')
s = traverse_freefem(exp).simplify()
print("only stiffener horizontal:",rep_derivatives(str(s.subs(sp.Function('Vy')(x, y),0).subs(sp.Function('Vx')(x, y),sp.Function('Vx')(x)).doit())
                                    ,lag.fields).replace(f'Derivative(Vx({x}), x)', 'dx(Vx)').replace(f'Derivative(Vx({x}), (x, 2))', 'dxx(Vx)')) # for vertical stiffener