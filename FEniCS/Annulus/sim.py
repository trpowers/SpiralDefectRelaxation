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


# %%
pi=np.pi
R=1/np.exp(1) #Radius of inclusion
R2 = 1


lsq = 0.001 #(K/A) in domain units



num_steps=90000;
N = 60;
frame_interval = 100;

alpha0_set = np.linspace(0.6, 1, 5)*np.pi
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


domain1 = domain2 = Circle(Point(0., 0.), R2)
domain2 = Circle(Point(0., 0.), R)
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


for phi_value in np.array([0.51*pi]):
    for j in range(m, m+1):
        # High-level FEniCS string expression using standard C++ logic
        cpp_expression = """
        ({
            const double pi_val = 3.14159265358979323846;
            std::complex<double> z(x[0], x[1]);
            double r2 = x[0]*x[0] + x[1]*x[1];
            double result = 0.0;

            if (r2 >= 1e-15) {
                std::complex<double> zp(r * std::cos(phi),  r * std::sin(phi));
                std::complex<double> zm(r * std::cos(phi), -r * std::sin(phi));

                double l = std::log(1.0 / rho);
                std::complex<double> exp_pow(0.0, -pi_val / l);

                std::complex<double> xi  = std::pow(z, exp_pow);
                std::complex<double> xip = std::pow(zp, exp_pow);
                std::complex<double> xim = std::pow(zm, exp_pow);

                std::complex<double> term1 = (xi - xip) / (xi - xim);
                std::complex<double> term2 = (xi - std::conj(xip)) / (xi - std::conj(xim));
                
                double arg_w = std::arg(term1 * term2);
                double log_z_conj_z = std::log(r2);

                result = 0.5 * arg_w + (alpha0 / 2.0) * log_z_conj_z / std::log(rho);
            }
            result; // Evaluates out of the structural block
        })
        """

        # FEniCS parses these keyword arguments automatically into C++ variables!
        alpha0 = alpha0_set[j]
        r = R**(0.5)
        theta = fe.Expression(cpp_expression,degree=2,rho=R,alpha0=alpha0, r=r, phi=phi_value)

        xdmf_file = XDMFFile(full_mesh.mpi_comm(),"./solutions/defects_solution_alpha0{0:.2f}_phi{1:.2f}_r{2:.2f}_KA{3:.4f}_rho{4:.2f}.xdmf".format(alpha0, phi_value, r, lsq, R))
        xdmf_file.parameters["flush_output"] = True
        xdmf_file.parameters["functions_share_mesh"] = True # Important!
        xdmf_file.parameters["rewrite_function_mesh"] = False # Important for static mesh!
        k_val = j*pi/(np.log(R/R2))
        m_val = np.log(R2)


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
            return on_boundary and sqrt(x[0]**2+x[1]**2)<R+tolbc #x.distance(p) < R+tol
        
        def outer_boundary(x, on_boundary):
            return on_boundary and sqrt(x[0]**2+x[1]**2)>R+tolbc #x.distance(p) > R+tol


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
