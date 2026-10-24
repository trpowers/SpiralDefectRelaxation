# %%
from __future__ import print_function
from fenics import *
from mshr import *
import matplotlib.pyplot as plt
from dolfin import *
from math import pi
import numpy as np
import random
import math
import ufl
import matplotlib
import sys
import os
import fenics as fe
import time

# %%
pi=np.pi
u1 = 1.5
u2 = 2
R2 = np.cosh(u2)
l_val = u2-u1



lsq = 0.001 #sqrt(K/A) in domain units



num_steps=30000;
N = 60;
frame_interval = 100;

alpha0_set = np.array([0.8, 0.9, 1])*np.pi
# %%
def get_eig(hes):
    mesh = hes.function_space().mesh()
    [eigL,eigR] = np.linalg.eig(hes.vector().get_local().reshape([mesh.num_vertices(),2,2]))
    ind = np.argsort(eigL,axis=1)[:,::-1]
    #print(eigR)
    eigL = np.take_along_axis(eigL,ind,axis=1)
    eigR = np.array([eigR[i][:,ind[i]] for i in range(len(eigR))])
    eigg = Function(TensorFunctionSpace(mesh,'CG',1,(2,2)))
    eig = Function(VectorFunctionSpace(mesh,'CG',1))
    eig.vector().set_local(eigL.flatten())
    eigg.vector().set_local(eigR.flatten())
    return eig, eigg

# %% [markdown]
# 

domain1 = Ellipse(Point(0., 0.), 1, np.sinh(u2)/np.cosh(u2))
domain2 = Ellipse(Point(0., 0.), np.cosh(u1)/np.cosh(u2), np.sinh(u1)/np.cosh(u2))
domain = domain1-domain2
full_mesh = generate_mesh(domain, N)
hmin=full_mesh.hmin();print(hmin)
hmax=full_mesh.hmax();print(hmax)
tol=hmin
tolbc=hmin

SF_full_mesh = FunctionSpace(full_mesh, 'CG', 1)
VS_full_mesh = VectorFunctionSpace(full_mesh, 'CG', 1)
#define function space
TS_full_mesh = TensorFunctionSpace(full_mesh, 'CG', 1,(2,2))
VE = VectorElement("Lagrange", full_mesh.ufl_cell(), 2)
FE = FiniteElement("Lagrange", full_mesh.ufl_cell(), 1)
W = FunctionSpace(full_mesh,MixedElement([FE,FE]))
# Define trial and test functions
w = Function(W)
dw = TrialFunction(W)
(Qxx,Qxy) = split(w)
(qxx,qxy) = TestFunctions(W)

m=int(sys.argv[1])


for phi_value in np.array([pi/6]):
    for j in range(m, m+1):
        # High-level FEniCS string expression using standard C++ logic
        cpp_expression = """
({
    const double pi_val = 3.14159265358979323846;
    double result = 0.0;
    
    // --- 1. Field Point: Calculate u, v for (x, y) ---
    double d1 = std::sqrt((x[0] + a)*(x[0] + a) + x[1]*x[1]);
    double d2 = std::sqrt((x[0] - a)*(x[0] - a) + x[1]*x[1]);
    
    double arg_acosh = (d1 + d2) / (2.0 * a);
    if (arg_acosh < 1.0) arg_acosh = 1.0; 
    double u_val = std::acosh(arg_acosh);
    
    double arg_acos = (d1 - d2) / (2.0 * a);
    if (arg_acos > 1.0) arg_acos = 1.0;
    if (arg_acos < -1.0) arg_acos = -1.0;
    
    double sign_y = (x[1] >= 0.0) ? 1.0 : -1.0;
    double v_val = sign_y * std::acos(arg_acos);
    
    // --- 2. Conformal mapping for field point (xi) ---
    // Note: Mathematica script changed (u - u1) to (u - u2)
    std::complex<double> I_comp(0.0, 1.0); 
    std::complex<double> exp_pow(0.0, -pi_val / l);
    
    std::complex<double> w_xi(u_val - u2, v_val);
    std::complex<double> xi = std::exp(exp_pow * w_xi);
    
    // --- 3. Defect Point: Calculate u_d, v_d for (r*cos(phi), r*sin(phi)) ---
    double xd = r * std::cos(phi);
    double yd = r * std::sin(phi);
    
    double d1_d = std::sqrt((xd + a)*(xd + a) + yd*yd);
    double d2_d = std::sqrt((xd - a)*(xd - a) + yd*yd);
    
    double arg_acosh_d = (d1_d + d2_d) / (2.0 * a);
    if (arg_acosh_d < 1.0) arg_acosh_d = 1.0;
    double u_d = std::acosh(arg_acosh_d);
    
    double arg_acos_d = (d1_d - d2_d) / (2.0 * a);
    if (arg_acos_d > 1.0) arg_acos_d = 1.0;
    if (arg_acos_d < -1.0) arg_acos_d = -1.0;
    
    double sign_yd = (yd >= 0.0) ? 1.0 : -1.0;
    double v_d = sign_yd * std::acos(arg_acos_d);

    // --- 4. Conformal mapping for defect (xiplus and ximinus) ---
    // Since ximinus evaluates (-yd), u_d remains the same but v_d becomes -v_d
    std::complex<double> wp(u_d - u2, v_d);
    std::complex<double> wm(u_d - u2, -v_d);
    
    std::complex<double> xiplus = std::exp(exp_pow * wp);
    std::complex<double> ximinus = std::exp(exp_pow * wm);
    
    // --- 5. Calculate defect field omegad -> thetad ---
    std::complex<double> term1 = (1.0 - xiplus/xi) / (1.0 - ximinus/xi);
    std::complex<double> term2 = (1.0 - std::conj(xiplus)/xi) / (1.0 - std::conj(ximinus)/xi);
    
    std::complex<double> omegad = (-I_comp / 2.0) * std::log(term1 * term2);
    double thetad = std::real(omegad);
    
    // --- 6. Calculate background field thetasp (Unchanged) ---
    double thetasp = alpha * (u2 - u_val) / (u2 - u1);
    
    // --- 7. Superposition ---
    result = thetasp + thetad;
    
    result;
})
"""
        alpha0 = alpha0_set[j]
        r=np.sqrt(np.cosh(u1)/np.cosh(u2))
# FEniCS Expression call passing the necessary Python variables to C++
        theta = fe.Expression(
            cpp_expression, 
            degree=2, 
            a=1/np.cosh(u2),# Focal length of the ellipses
            u1=u1,          # Inner ellipse coordinate parameter
            u2=u2,          # Outer ellipse coordinate parameter
            l=l_val,        # Length scale (u2 - u1)
            r=np.sqrt(np.cosh(u1)/np.cosh(u2)),        # Radial position of the defect in the reference plane
            phi=phi_value,  # Angular position of the defect
            alpha=alpha0 # Boundary condition twist magnitude
        )


        xdmf_file = XDMFFile(full_mesh.mpi_comm(),"./solutions/ellipse_solution_alpha0{0:.2f}_phi{1:.2f}_r{2:.2f}_KA{3:.4f}_u1{4:.2f}_u2{5:.2f}.xdmf".format(alpha0, phi_value, r, lsq, u1, u2))
        xdmf_file.parameters["flush_output"] = True
        xdmf_file.parameters["functions_share_mesh"] = True # Important!
        xdmf_file.parameters["rewrite_function_mesh"] = False # Important for static mesh!


        Qxxi = project(Expression("cos(2*w)", w=theta, degree =1), SF_full_mesh)
        Qxyi = project(Expression("sin(2*w)", w=theta, degree =1), SF_full_mesh)

        Qi = as_tensor([[Qxxi, Qxyi], [Qxyi, -Qxxi]])

        Qi_projected = project(Qi, TensorFunctionSpace(full_mesh, 'CG', 1, (2,2)))


        Qn = project(as_tensor([[Qxxi, Qxyi], [Qxyi, -Qxxi]]), TS_full_mesh)
        dt =0.01

        t = 0

        Q = as_tensor([[Qxx, Qxy], [Qxy, -Qxx]])

        q = as_tensor([[qxx, qxy], [qxy, -qxx]])  


        F1 =(1-dt)*inner(Q,q)*dx + lsq*dt*inner(grad(q), grad(Q))*dx -inner(Qn,q)*dx + dt*inner(dot(dot(Q,Q), Q),q)*dx  #Q equation 


        F=F1

        J = derivative(F, w, dw)

        def inner_boundary(x, on_boundary):
            return on_boundary and sqrt(x[0]**2+x[1]**2)<=np.sqrt(np.cosh(u1)/np.cosh(u2)) #x.distance(p) < R+tol
        
        def outer_boundary(x, on_boundary):
            return on_boundary and sqrt(x[0]**2+x[1]**2)>np.sqrt(np.cosh(u1)/np.cosh(u2)) #x.distance(p) > R+tol


        Qxx_inner= DirichletBC(W.sub(0), Constant(np.cos(2*(alpha0))), inner_boundary)
        Qxy_inner= DirichletBC(W.sub(1), Constant(np.sin(2*(alpha0))), inner_boundary)

        Qxx_outer= DirichletBC(W.sub(0), Constant(1), outer_boundary)
        Qxy_outer= DirichletBC(W.sub(1), Constant(0), outer_boundary)




        bc = [Qxx_inner, Qxy_inner, Qxx_outer, Qxy_outer] #Qxx_wall, Qxy_wall removed
        S_zero = np.empty((0, 2))

        eg, ev =get_eig(project(Qn))
        nn = as_vector((ev[0, 0], ev[1, 0]))
        S=sqrt(Qxyi**2+Qxxi**2)
        S_func = project(S, SF_full_mesh)
        nn_func = project(nn, VS_full_mesh)
        Qxx_func = project(Qxxi, SF_full_mesh)
        Qxy_func = project(Qxyi, SF_full_mesh)
        S_func.rename("S_order_param", "S") # Renames the field to "S_order_param"
        nn_func.rename("director", "nn") # Renames the field to "director"
        Qxx_func.rename("Qxx", "qxx")
        Qxy_func.rename("Qxy", "qxy")
        xdmf_file.write(S_func, t)
        xdmf_file.write(nn_func, t)
        xdmf_file.write(Qxx_func, t)
        xdmf_file.write(Qxy_func, t)
        time_start = time.time()
    # %%
        for n in range(num_steps):
            # Solve variational problem
            w_guess_for_this_step = Function(w.function_space())
            w_guess_for_this_step.assign(w)
            problem = NonlinearVariationalProblem(F, w, bc, J)
            solver  = NonlinearVariationalSolver(problem)
            prm = solver.parameters
            if n==0:
                default_solver = prm['newton_solver']['linear_solver']
            rel = 0.9
            nmi=100
            abs_tol = 1E-8
            convergence_flag = 0
            # Newton Solver Parameters 
            while(convergence_flag==0):
                w.assign(w_guess_for_this_step)
                prm["newton_solver"]["maximum_iterations"] = nmi
                prm["newton_solver"]["relaxation_parameter"] = rel
                # prm['newton_solver']['relative_tolerance'] = 1E-6
                prm['newton_solver']['absolute_tolerance'] = abs_tol

                try:
                    num_iterations, converged = solver.solve()
                except RuntimeError as e:
                    print(f"  RuntimeError during solve: {e}")
                    converged = False
                
                if converged:
                    Qxx,Qxy = split(w)

                    # # Update current time
                    t = round(t+dt,4);print('t='+str(t))

                    eg, ev =get_eig(project(Q))
                    nn = as_vector((ev[0, 0], ev[1, 0]))
                    S=sqrt(Qxy**2+Qxx**2)

                    Q_current_ufl = as_tensor([[Qxx, Qxy], [Qxy, -Qxx]])
                    Qn.assign(project(Q_current_ufl, TS_full_mesh))        
                    #Qn=Q
                    sys.stdout.flush()
                    if (n%frame_interval==0):
                        S_func = project(S, SF_full_mesh)
                        nn_func = project(nn, VS_full_mesh)
                        Qxx_func = project(Qxx, SF_full_mesh)
                        Qxy_func = project(Qxy, SF_full_mesh)
                        S_func.rename("S_order_param", "S") # Renames the field to "S_order_param"
                        nn_func.rename("director", "nn") # Renames the field to "director"
                        Qxx_func.rename("Qxx", "qxx")
                        Qxy_func.rename("Qxy", "qxy")
                        xdmf_file.write(S_func, t)
                        xdmf_file.write(nn_func, t)
                        xdmf_file.write(Qxx_func, t)
                        xdmf_file.write(Qxy_func, t)
                    convergence_flag = 1
                    prm['newton_solver']['linear_solver']=default_solver
                else:
                    prm['newton_solver']['linear_solver']='mumps'
                    abs_tol = 1E-5
                    if rel<0.2:
                        rel = rel*0.95
                        nmi += 100
                    else:
                        rel = rel*0.9
        print(time.time()-time_start, num_steps)
