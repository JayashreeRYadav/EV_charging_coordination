using JuMP
using Ipopt
using Gurobi
using Plots
using LaTeXStrings
using MAT


#=
M = 100;  # Number of iterations 
objVal = zeros(M,1);

objVal = zeros(M,1);
gamma_val = zeros(M,1);
for ii=1:M
    gamma_val[ii] = (ii-1)*2/M;
end
xval = zeros(M,1);
#yval = zeros(M,1);
#uval = zeros(M,1);

eta = 0.9;
K = 2; # No of PVs connected
P_max = 50; #kW
t = 1:1:24;
D = [24, 22, 21, 20, 21, 26, 31, 33, 35, 36, 37, 36, 35, 35.5, 34.5, 34, 34.7, 35, 36, 35, 34, 32, 30, 25];
p1 = zeros(24,M);
p2 = zeros(24,M)
vt = zeros(M,1); 
vt1 = 0;
alpha = 0.5;
for m in 1:M
    for t in 1:24
        vt[m] = (D[t] + (p1[t,m]+p2[t,m]))^2 + gamma_val[m]*((p1[t,m]+p2[t,m])-eta*P_max);
    
        p1[t,m]= p1[t,m] -alpha*vt[m];
        p2[t,m] = p2[t,m]-alpha*vt[m];
    end
end

m = Model(Gurobi.Optimizer);
p1t1 = @variable(m);
p1t2 = @variable(m);
p2t1 = @variable(m);
p2t2 = @variable(m);
#y = @variable(m);
d1 = 14; d2 = 12
@constraint(m, p1t1+p1t2 == 10);
@constraint(m, p2t1+p2t2 == 10);
@constraint(m, p1t1>=0);
@constraint(m, p1t2>=0);
@constraint(m, p2t1>=0);
@constraint(m, p2t2>=0);
@constraint(m, d1+d2+p1t1+p1t2+p2t2+p2t1 <= 90);
@objective(m, Min, (d1+d2+p1t1+p1t2+p2t2+p2t1)^2);

optimize!(m)
ObjValue = objective_value(m);
p1t1Value=value(p1t1);
p1t2Value=value(p1t2);
p2t1Value=value(p2t1);
p2t2Value=value(p2t2);
println("C(x) = ",ObjValue);
println("xPrimal = ",p1t1Value);
println("xPrimal = ",p1t2Value);
println("xPrimal = ",p2t1Value);
println("xPrimal = ",p2t2Value);
=#
#=
#for 1 EV
N = 10; # time steps
eta = 0.9; # feeder efficiency
Pmax = 50; # max feeder capacity
pub = 2;
plb = 0;
#T = 10;
D = [36, 24, 22, 21, 20, 21, 26, 31, 33, 35]; # base load
m = Model(Gurobi.Optimizer);#direct_model(Gurobi.Optimizer());
p=@variable(m,[1:N]);
@objective(m, Min, sum((D[n] + p[n])^2 for n = 1:N))
@constraint(m, p .<= pub)
@constraint(m, p .>= plb);
@constraint(m, cc, sum(p[n] for n = 1:N) == 10)
#@constraint(m, [D[n] + p[n] <= eta * Pmax for n in 1:N])
for n in 1:N
    @constraint(m, D[n] + p[n] <= eta * Pmax)
end
grb_model = JuMP.backend(m);
optimize!(m)
ObjValue = objective_value(m);
p_value = value.(p);
println("1 EV and 10 time steps");
println("C(x) = ",ObjValue);
println("xPrimal = ",p_value);
=#


# for K EVs
# for N EVs
N = 5; # No of EVs
T = 10; # time steps
eta = 0.9; # feeder efficiency
Pmax = 50; # max feeder capacity
pub = 2;
plb = 0;
Umax = 10;
#T = 10;
D = [36, 24, 22, 21, 20, 21, 26, 31, 33, 35]; # base load
solveFlag=1;
if(solveFlag==1)
    m = Model(Ipopt.Optimizer);#direct_model(Gurobi.Optimizer());
elseif(solveFlag==0)
    m = Model(Gurobi.Optimizer);#direct_model(Gurobi.Optimizer());
end
@variable(m, p[1:N, 1:T])
@objective(m, Min, sum((D[t] + sum(p[n, t] for n = 1:N))^2 for t = 1:T))
@constraint(m, p .<= pub)
@constraint(m, p .>= plb);
#@constraint(m, cc, (sum(p[n, t] for t = 1:T) == 10) for n = 1:N)
#@constraint(m, [D[t] + p[t] <= eta * Pmax for t in 1:T])
for n in 1:N
    @constraint(m, sum(p[n, t] for t = 1:T) == Umax);
end
for t = 1:T
    @constraint(m, D[t] + sum(p[n, t] for n = 1:N) <= eta * Pmax);
end
grb_model = JuMP.backend(m);
optimize!(m)
ObjValue = objective_value(m);
p_value = value.(p);
println("10 EV and 10 time steps");
println("C(x) = ",ObjValue);
println("xPrimal = ",p_value);
if solveFlag==0
    matwrite("Data_2norm_Gurobi.mat", Dict("ev_charge_Gurobi" => p_value, "D" => D));
elseif solveFlag == 1
    matwrite("Data_2norm_Ipopt.mat", Dict("ev_charge_ipopt" => p_value,  "D" => D));
end
plot(title="EV charging Distribution 2 norm", ylabel="Value", xlabel="time step")
plot!(1:T, D, label="Base load")
plot!(1:T, sum(p_value[n,:] for n = 1:N), label="EV charging load")
plot!(1:T, D + sum(p_value[n,:] for n = 1:N), label="Base load + EV")

#=
println("\nSolving dualDecomp")
K = 10
T = 10
eta = 0.9
Pmax = 50
pub = 2
plb = 0
D = [36, 24, 22, 21, 20, 21, 26, 31, 33, 35]
itr_max = 100;
epsTol = 1e-6;
λ = zeros(K, itr_max);  
μ = zeros(T, itr_max);  
global itr = 1;
# Dual decomposition loop
for itr = 1: itr_max-1 #eps[itr] > epsTol  
    #println("Iteration ", itr)
    
   
    #for t in 1:T
        m1 = Model(Gurobi.Optimizer);
        @variable(m1, pkt[1:K, 1:T]);
        @constraint(m1, pkt .>= plb)
        @constraint(m1, pkt .<= pub)
        @objective(m1, Min, sum((D[t] + sum(pkt[k, t] for k = 1:K))^2 for t = 1:T) + sum( μ[t, itr] * (+ D[t] + sum(pkt[k, t] for k = 1:K) - eta * Pmax)  for t = 1:T )+ sum( λ[k, itr] * (sum(pkt[k, :]) - 10) for k = 1:K )); 
        optimize!(m1);
        β = 1/(itr)^2;
        for t = 1:T
            μ[t, itr+1] = μ[t, itr] + β .* (D[t]+ sum(value.(pkt[:,t])) - eta * Pmax);
        end
         # obj_values = objective_value.(m1)
        #p_values = value.(pkt)
   # end
        α = 1/(itr)^2;
         for k = 1:K   
            λ[k, itr+1] = λ[k, itr] + α .* (sum(value.(pkt[k, :])) - 10)
         end
    #epsTol = λ[itr+1] - λ[itr];
    itr = itr + 1;
end

# Print results
#println("Objective value: ", objective_value(m_master))
#println("Primal variable values: ", value.(p))
=#

## Extracting suboptimal solution/Check if strong or weak minimum
K=10;ϵ = 0;
solution_matrices = []
ObjValue_subopt = zeros(K,1)
for k = 1:K
c = rand(N,T);
m2 = Model(Gurobi.Optimizer);
@variable(m2, p2[1:N, 1:T]);
@objective(m2, Min, sum(c[n,t] * p2[n, t] for n in 1:N, t in 1:T))
@constraint(m2, p2 .<= pub);
@constraint(m2, p2 .>= plb);
@constraint(m2, sum((D[t] + sum(p2[n, t] for n = 1:N))^2 for t = 1:T) <= ObjValue+ϵ);
for n in 1:N
    @constraint(m2, sum(p2[n, t] for t = 1:T) == Umax);
end
for t = 1:T
    @constraint(m2, D[t] + sum(p2[n, t] for n = 1:N) <= eta * Pmax);
end
optimize!(m2);
ObjValue_subopt[k] = objective_value(m2);
push!(solution_matrices, value.(p2));
# matwrite("solution_matrices_$k.mat", Dict("solution_matrix" => value.(p)))
end
filename = "solution_matrices.mat"
solution_dict = Dict("solution_matrix_$i" => solution_matrices[i] for i in 1:length(solution_matrices))
matwrite(filename, solution_dict)