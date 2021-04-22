function PXPExactDynamics
% id/dt \rho = [H_full, \rho]

%%parameters
global imagunit tolerance
imagunit = sqrt(-1);
tolerance = 1.0e-12;
lambda = 1;
total_iter = 1000;
dt = 0.01;
tem = 1; %temperature
beta = 1/tem;
N_sites_A = 2; %num of sites in subsystem A
N_sites_B = 8; %num of sites in subsystem B
N_sites = N_sites_A + N_sites_B;

%%bases in real space
bases_AB = Bases(N_sites);
bases_AB_constrained = Bases_constrained(bases_AB, N_sites);
disp(size(bases_AB_constrained, 1))
[bases_AB_constrained_1, bases_AB_constrained_2] = split_bases_AB(bases_AB_constrained, N_sites_A);
bases_AB_constrained_12 = [bases_AB_constrained_1 ; bases_AB_constrained_2];
dim_subspace1 = size(bases_AB_constrained_1, 1)
dim_subspace2 = size(bases_AB_constrained_2, 1)

%%Hamiltonian of the whole system, including A and B
H_AB_12_real = Hamiltonian_full(bases_AB_constrained_12, N_sites_A, lambda);
H_AB_1_real = Hamiltonian_full(bases_AB_constrained_1, N_sites_A, lambda);
H_AB_2_real = Hamiltonian_full(bases_AB_constrained_2, N_sites_A, lambda);
[eigenstates_1,eigen_energies_mat_1] = eig(H_AB_1_real);
[eigenstates_2,eigen_energies_mat_2] = eig(H_AB_2_real);
eigenstates_1_expanded = [eigenstates_1; zeros(dim_subspace2, dim_subspace1)];
eigenstates_2_expanded = [zeros(dim_subspace1, dim_subspace2); eigenstates_2];
bases_energy_12 = [eigenstates_1_expanded, eigenstates_2_expanded];
H_AB_12_energy = bases_energy_12' * H_AB_12_real * bases_energy_12;

[eigenstates_12,eigen_energies_mat_12] = eig(H_AB_12_real);
rho_full = zeros(size(H_AB_12_real));
for i=1:size(H_AB_12_real,1)
    rho_full_thermal(i,i) = exp(-beta*eigen_energies_mat_12(i,i));
end
bases_energy_12' * rho_full_thermal * bases_energy_12


%%initial density matrix \rho_{full}
rho11_energy = zeros(size(eigen_energies_mat_1,1));
for i=1:size(eigen_energies_mat_1,1)
    rho11_energy(i,i) = exp(-beta * eigen_energies_mat_1(i,i));
end
rho22_energy = zeros(size(eigen_energies_mat_2,1));
for i=1:size(eigen_energies_mat_2,1)
    rho22_energy(i,i) = exp(-beta * eigen_energies_mat_2(i,i));
end
rho_full_energy = blkdiag(rho11_energy, rho22_energy);
rho_full_energy = rho_full_energy./trace(rho_full_energy);

%%Exact dynamics
entropy_full_arr = [];
entropy_full_arr = [entropy_full_arr, Entropy(rho_full_energy)];
rho12block_11_arr = [];
rho12block = rho_full_energy(1:dim_subspace1, dim_subspace1+1:end);
rho12block_11_arr = [rho12block_11_arr, rho12block(1,1)];
rho11block_12_arr = [];
rho11block_12_arr = [rho11block_12_arr, rho_full_energy(1,2)];
params = {H_AB_12_energy};
for n_iter = 2 : total_iter
    rho_full_energy = RungeKutta_oneMatrix(rho_full_energy, @rhofullEvolution, params, dt);
    rho11block_12_arr = [rho11block_12_arr, rho_full_energy(1,2)];
    entropy_full_arr = [entropy_full_arr, Entropy(rho_full_energy)];
    rho12block = rho_full_energy(1:dim_subspace1, dim_subspace1+1:end);
    rho12block_11_arr = [rho12block_11_arr, rho12block(1,1)];
end

rho_full_energy

%%plots
tvector=0:dt:(total_iter-1)*dt;
fig_num = 1;

%entropy_full vs time
figure(fig_num);
fig_num = fig_num + 1;
plot(tvector, real(entropy_full_arr), 'blue--');
%hold on
%plot(tvector, imag(entropy_full_arr), 'red-.')
title("Entropy vs time")

%rho elements vs time
figure(fig_num);
fig_num = fig_num + 1;
plot(tvector, abs(rho12block_11_arr), 'b--');
title("rho12(1,1) magnitude vs time");
figure(fig_num);
fig_num = fig_num + 1;
plot(tvector, Complex_angle(rho12block_11_arr), 'r-.');
title("rho12(1,1) angle vs time")

figure(fig_num);
fig_num = fig_num + 1;
plot(tvector, abs(rho11block_12_arr), 'b--');
title("rho11(1,2) magnitude vs time");
figure(fig_num);
fig_num = fig_num + 1;
plot(tvector, Complex_angle(rho11block_12_arr), 'r-.');
title("rho11(1,2) angle vs time")


end

function complex_angles = Complex_angle(z)
%complex angles of a complec vector
complex_angles = zeros(size(z));
for i = 1:size(z,2)
    complex_angles(i) = atan(imag(z(i))/real(z(i)));
end
end
  


%%bases
function [bases_B_1, bases_B_2]=split_bases_B(bases_B)
%split the bases of B into two subspaces with the first spin being up or
%down
bases_B_1 = [];
bases_B_2 = [];
n_states = size(bases_B, 1);
for i=1:n_states
    basis_i = bases_B(i,:);
    if basis_i(1) == -1
        bases_B_1 = [bases_B_1;basis_i];
    else
        bases_B_2 = [bases_B_2;basis_i];
    end
end
end

function bases_a = Bases(N_sites)
%to find the bases in subsystem A or subsystem B
%spin down = -1
%spin up = 1
sites_vec = 1:N_sites;
total = 2^(N_sites);
bases_a = -ones(total,N_sites);
start = 1;
for i=1:N_sites
    spinup_ind = nchoosek(sites_vec,i);
    size_spinup_ind = size(spinup_ind,1);
    for j=1:size_spinup_ind
        bases_a(start+j,spinup_ind(j,:))=1;
    end
    start = start + size_spinup_ind;
end
end

function bases_a_constrained = Bases_constrained(bases_a, N_sites)
%filter bases_a to the constrained space, which doesn't allow two adjacent
%up spins.
dim = size(bases_a,1);
bases_a_constrained = [];
for i=1:dim
    basis = bases_a(i,:);
    flag = true;
    for j=1:N_sites-1
        if basis(j)>0 && basis(j+1)>0
            flag = false;
            break;
        end
    end
    if flag 
        bases_a_constrained = [bases_a_constrained;basis];
    end
end
end

function [bases_AB_1, bases_AB_2] = split_bases_AB(bases_AB, N_sites_A)
%splite the bases of the whole system, including A and B, into subspace 1 and 2. 
%In subpace 1, spin on site N_sites_A + 1 is down, 
%while in subspace 2, spin on site N_sites_A + 1 is up. 
%downspin : -1
%upspin: +1
bases_AB_1 = [];
bases_AB_2 = [];
n_states = size(bases_AB, 1);
for i=1:n_states
    basis_i = bases_AB(i,:);
    if basis_i(N_sites_A + 1) == -1
        bases_AB_1 = [bases_AB_1;basis_i];
    else
        bases_AB_2 = [bases_AB_2;basis_i];
    end
end
end

%%Hamiltonian
function H_b = Hamiltonian_b(bases)
dimH = size(bases,1);
H_b = zeros(dimH,dimH);
for i=1:dimH
    for j=i+1:dimH
        basis_i = bases(i,:);
        basis_j = bases(j,:);
        H_b(j,i) = PXP_b(basis_i, basis_j);
    end
end
H_b = H_b + H_b';
end

function Hji = PXP_b(basis_i, basis_j)
%Hamiltonian of subsystem B
%H = P_1 X_2 P_3 + ... + P_{L-1} X_L
%used by Hamiltonian_b
N_sites = size(basis_i, 2);
diff = basis_i - basis_j;
diff_elements = find(diff);

Hji = 0;
if size(diff_elements,2)==1 
    if diff_elements(1) == N_sites
        if basis_i(N_sites-1) == -1
            Hji = 1;
        end
    else       
        xi = diff_elements(1);
        if basis_i(xi-1)==-1 && basis_i(xi+1)==-1
            Hji = 1;
        end
    end
end
end

function H_a = Hamiltonian_a(bases)
% to find the matrix form of Hamiltonian_a in the bases
% functions used: PXP
dimH = size(bases,1);
H_a = zeros(dimH,dimH);
for i=1:dimH
    for j=i+1:dimH
        basis_i = bases(i,:);
        basis_j = bases(j,:);
        H_a(j,i) = PXP_a(basis_i,basis_j);
    end
end
H_a = H_a + H_a';
end

function Hji = PXP_a(basis_i,basis_j)
% the matrix element between basis_i and basis_j
% H = X_1 P_2 + P_1 X_2 P_3 + ... + P_{L-1} X_L
% where L=N_sites
% used by Hamiltonian_a
N_sites = size(basis_i, 2);
diff = basis_i - basis_j;
diff_elements = find(diff);

Hji = 0;
if size(diff_elements,2)==1 
    if diff_elements(1)==1
        if basis_i(2)==-1
            Hji = 1;
        end
    elseif diff_elements(1) == N_sites
        if basis_i(N_sites-1) == -1
            Hji = 1;
        end
    else       
        xi = diff_elements(1);
        if basis_i(xi-1)==-1 && basis_i(xi+1)==-1
            Hji = 1;
        end
    end
end
end

function H_full = Hamiltonian_full(bases, N_sites_A, lambda)
% to find the matrix form of Hamiltonian_full in the bases
% functions used: PXP_full
dimH = size(bases,1);
H_full = zeros(dimH,dimH);
for i=1:dimH
    for j=i+1:dimH
        basis_i = bases(i,:);
        basis_j = bases(j,:);
        H_full(j,i) = PXP_full(basis_i,basis_j, N_sites_A, lambda);
    end
end
H_full = H_full + H_full';
end

function Hji = PXP_full(basis_i,basis_j, N_sites_A, lambda)
% the matrix element between basis_i and basis_j
% H = X_1 P_2 + P_1 X_2 P_3 + ...+ \lambda*P_L X_{L+1} P_{L+2}  ... + P_{L-1} X_L
% where L=N_sites
% used by Hamiltonian_a
Lp1 = N_sites_A + 1;
N_sites = size(basis_i, 2);
diff = basis_i - basis_j;
diff_elements = find(diff);

Hji = 0;
if size(diff_elements,2)==1 
    if diff_elements(1)==1
        if basis_i(2)==-1
            Hji = 1;
        end
    elseif diff_elements(1) == N_sites
        if basis_i(N_sites-1) == -1
            Hji = 1;
        end
    elseif diff_elements(1) == Lp1
        if basis_i(Lp1-1)==-1 && basis_i(Lp1+1)==-1
            Hji = lambda;
        end
    else       
        xi = diff_elements(1);
        if basis_i(xi-1)==-1 && basis_i(xi+1)==-1
            Hji = 1;
        end
    end
end
end



%%write operators in the form of matrices
function PL_matrix = PL(from_subspace,to_subspace)
%matrix form of P_L, going from one subspace to the other subspace of
%subsystem A
%either from_subpace or to_subspace must be in complate shape, i.e. the
%states in subspace 2 of subsystem A cannot be truncated (the last spin)
n_sites = size(from_subspace, 2);
dim_from = size(from_subspace, 1);
dim_to = size(to_subspace, 1);
PL_matrix = zeros(dim_to, dim_from);
for fromi = 1:dim_from
    state_from = from_subspace(fromi,:);
    for toj = 1:dim_to
        state_to = to_subspace(toj,:);
        if state_from(n_sites)==-1 && state_to(n_sites)==-1
            state_abs_diff = abs(state_from - state_to);
            if sum(state_abs_diff) == 0
                PL_matrix(toj,fromi) = 1;
            end
        end
    end
end
end

function XLp1PLp2_matrix = XLp1PLp2(from_subspace, to_subspace)
%matrix form of X_{L+1}P_{L+2}, going from one subspace to the other subspace of
%subsystem A
dim_from = size(from_subspace, 1);
dim_to = size(to_subspace, 1);
XLp1PLp2_matrix = zeros(dim_to, dim_from);
for fromi = 1:dim_from
    state_from = from_subspace(fromi,:);
    for toj = 1:dim_to
        state_to = to_subspace(toj,:);
        if state_from(2) == -1 && state_to(2) == -1
            state_abs_diff = abs(state_from - state_to);
            if state_abs_diff(1) == 2 && sum(state_abs_diff) == 2
                XLp1PLp2_matrix(toj,fromi) = 1;
            end
        end
    end
end
end

%%Runge-Kutta, evolution of density matrices
function rho_next = RungeKutta_oneMatrix(rho, rho_der_fun, params, dt)
%4th order Runge-Kutta
%one matrix 
k1 = rho_der_fun(rho, params);
k2 = rho_der_fun(rho + 0.5*dt*k1, params);
k3 = rho_der_fun(rho + 0.5*dt*k2, params);
k4 = rho_der_fun(rho + dt*k3, params);
rho_next = rho + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

function rhofull_der = rhofullEvolution(rho_full, params)
%evolution of rho_full
%params is a cell, {H_full}
global imagunit
H_full = params{1};
rhofull_der = -imagunit * commutator(H_full, rho_full);
end

function rho11a_der = rho11aEvolution0(rho11a, params)
%evolution of rho11_a
%params  is a cell, {H_a_1}
global imagunit
H_a_1 = params{1};
rho11a_der = -imagunit * commutator(H_a_1, rho11a);
end

function rho11full_der = rho11fullEvolution0(rho11full, params)
%evolution of full rho11 with lambda=0
%params is a cell, {H_a_1_full, H_b_1_full}
global imagunit
H_a_1_full = params{1};
H_b_1_full = params{2};
rho11full_der = -imagunit * commutator(H_a_1_full + H_b_1_full, rho11full);
end

function rho22full_der = rho22fullEvolution0(rho22full, params)
%evolution of full rho22 with lambda=0
%params is cell, {H_a_2_full, H_b_2_full}
global imagunit
H_a_2_full = params{1};
H_b_2_full = params{2};
rho22full_der = -imagunit * commutator(H_a_2_full + H_b_2_full, rho22full);
end

function rho12full_der = rho12fullEvolution0(rho12full, params)
%evolutioin of rho12 with lambda=0
%params is a cell, {H_a_1_full, H_a_2_full, H_b_1_full, H_b_2_full}
global imagunit
H_a_1_full = params{1};
H_a_2_full = params{2};
H_b_1_full = params{3};
H_b_2_full = params{4};
rho12full_der = -imagunit * (H_a_1_full * rho12full - rho12full * H_a_2_full + H_b_1_full * rho12full - rho12full * H_b_2_full);
end

function rho21full_der = rho21fullEvolution0(rho21full, params)
%evolutioin of rho12 with lambda=0
%params is a cell, {H_a_1_full, H_a_2_full, H_b_1_full, H_b_2_full}
global imagunit
H_a_1_full = params{1};
H_a_2_full = params{2};
H_b_1_full = params{3};
H_b_2_full = params{4};
rho21full_der = -imagunit * (H_a_2_full * rho21full - rho21full * H_a_1_full + H_b_2_full * rho21full - rho21full * H_b_1_full);
end


%%operation on matrices
function com = commutator(A, B)
com = A * B - B * A;
end

function anticom = anticommutator(A, B)
anticom = A * B + B * A;
end

function matrix_B = traceOverA(matrix_full, dim_A, dim_B)
%trace over B of the tensor product of two square matrices
%trace over A: sum of all diagonal dim_B-by-dim_B blocks in the full
%matrix
matrix_B = zeros(dim_B);
for n=0:dim_A-1
    range = n*dim_B+1 : (n+1)*dim_B;
    matrix_B  = matrix_B + matrix_full(range,range);
end
end

function matrix_A = traceOverB(matrix_full, dim_A, dim_B)
%trace over A of the tensor product of two square matrices
%trace over B: (i,j) element of matrix_A is equal to the trace of the
%(i,j) block of dimension dim_B-by-dim_B
matrix_A = zeros(dim_A);
for i=1:dim_A
    ri = (i-1)*dim_B+1 : i*dim_B;
    for j=1:dim_A
        rj = (j-1)*dim_B+1 : j*dim_B;
        matrix_A(i,j) = trace(matrix_full(ri,rj));
    end
end
end

function entropy = Entropy(densityMatrix)
% Boltzmann constant k_B = 1
[eigenstates,diagD] = eig(densityMatrix);
diagd = diag(diagD); 
dim = size(diagD,1);
entropy = 0;
for i = 1:dim
    entropy = entropy - diagd(i).*log(diagd(i));
end
end

