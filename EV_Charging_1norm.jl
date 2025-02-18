using JuMP
using Ipopt
using Gurobi
using Plots
using MAT
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
@objective(m, Min, sum((D[t] + sum(p[n, t] for n = 1:N)) for t = 1:T))
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
println("1 EV and 10 time steps");
println("C(x) = ",ObjValue);
println("xPrimal = ",p_value);
if solveFlag==0
    matwrite("Data_1norm_Gurobi.mat", Dict("ev_charge_Gurobi" => p_value, "D" => D));
elseif solveFlag == 1
    matwrite("Data_1norm_Ipopt.mat", Dict("ev_charge_ipopt" => p_value,  "D" => D));
end

plot(title="EV charging Distribution 1 norm", ylabel="Value", xlabel="time step")
plot!(1:T, D, label="Base load")
plot!(1:T, sum(p_value[n,:] for n = 1:N), label="EV charging load")
plot!(1:T, D + sum(p_value[n,:] for n = 1:N), label="Base load + EV")

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
    @constraint(m2, sum((D[t] + sum(p2[n, t] for n = 1:N)) for t = 1:T) <= ObjValue+ϵ);
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
filename = "solution_matrices_1norm.mat"
solution_dict = Dict("solution_matrix_$i" => solution_matrices[i] for i in 1:length(solution_matrices))
matwrite(filename, solution_dict)