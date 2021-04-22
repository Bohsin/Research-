function PXPLindblad
% Lindblad dynamics of PXP model 
% subsystem A:bases_a, hamiltonian_a, A_1 = P_{L-1}X_{L}, A_2 = P_{L}  
% bath:correlation functions
% parameters:
% N_sites -- length L of subsystem A, or the total number of sites in A

clear all
%%%parameters
global P_lp1_therm tol_energy T gamma imagunit;
imagunit=sqrt(-1);
P_lp1_therm = 0.3;
tol_energy = 1e-9;
T=0.5;
gamma = 10;
N_sites = 4;

dt = 0.15;
iteration_max = 8000;

bases_a = Bases_a(N_sites)
bases_a_constrained = Bases_a_constrained(bases_a, N_sites)
%H_a = Hamiltonian_a(bases_a,N_sites);
H_a_prime = Hamiltonian_a_prime(bases_a_constrained,N_sites)

dim = size(H_a_prime,1)

%A1 = zeros(dim,dim);
%for i=1:dim
%    for j=1:dim
%        A1(j,i) = A1_between_real_space_bases(bases_a,i,j,N_sites);
%    end
%end
%A1
%nonzeros(A1 ~= A1')

%%% two sectors of H_a
% [bases_a_lup, bases_a_ldown] = Symmetry_Zl(bases_a,N_sites);
% bases_a_lup
% H_a_lup = Hamiltonian_a(bases_a_lup,N_sites);
% bases_a_ldown
% H_a_ldown = Hamiltonian_a(bases_a_ldown,N_sites);
% 
% [eigenstates,D] = eig(H_a);
% energies = diag(D)'
% %eigenstates
% [eigenstates_lup,D_lup] = eig(H_a_lup);
% energies_lup = diag(D_lup)'
% %V_lup
% [eigenstates_ldown,D_ldown] = eig(H_a_ldown);
% energies_ldown = diag(D_ldown)'
% %V_ldown


[eigenstates_prime,H0matrix] = eig(H_a_prime);
energies_prime = diag(H0matrix)';

% A1 = zeros(dim,dim);
% for i=1:dim
%    for j=i:dim
%        %disp(['i=',num2str(i),'j=',num2str(j)])
%        A1(j,i) = A1_between_eigenstates(eigenstates,bases_a,i,j,N_sites);
%        A1(i,j) = conj(A1(j,i));
%        %pause
%    end
% end
% A1
% %diag(A1)
% %nonzeros(A1 ~= A1')


% figure(1)
% h1=histogram(energies_lup);
% h1.BinWidth = min(nonzeros(diff(energies_lup)))/5;
% xlabel("energy")
% ylabel("degeneracy")
% title(['Spin up on the site L',' (L=',num2str(N_sites),')'])
% 
% figure(2)
% h2=histogram(energies_ldown);
% h2.BinWidth = min(nonzeros(diff(energies_ldown)))/5;
% xlabel("energy")
% ylabel("degeneracy")
% title(['Spin down on the site L',' (L=',num2str(N_sites),')'])
% 
% figure(3)
% h3=histogram(energies);
% h3.BinWidth = min(nonzeros(diff(energies)))/5;
% xlabel("energy")
% ylabel("degeneracy")


figure(1)
h2=histogram(energies_prime);
h2.BinWidth = min(nonzeros(diff(energies_prime)))/5;
xlabel("energy")
ylabel("degeneracy")
title(['H_a','(prime)'])


% energy_deg = Degeneracy_spectra(energies_prime, dim)
% 
eigenstate_decomp = Decompose_eigenpairs(energies_prime, dim);
 
% size_ed = size(eigenstate_decomp.energies,2);
% for i=1:size_ed
%     eigenstate_decomp.energies(i)
%     eigenstate_decomp.eigen_pairs{i}
% end


%%%The Gamma matrix

%%%Matrices A: A1, A2
A2_matrix = A2_Matrix(eigenstates_prime,bases_a_constrained,N_sites);


%%%initial density matrix of subsystem A
%pure state
%mixed state
ini_ind = 8;
ini_state = zeros(dim,1);
ini_state(ini_ind) = 1;
ini_state_eigen_basis = eigenstates_prime'*ini_state;


rho = ini_state_eigen_basis*ini_state_eigen_basis';


energy_exp_vec = zeros(1,dim);
energy_exp_vec(1) = trace(rho*H0matrix);
entropy = zeros(1,dim);
entropy(1) = Entropy(rho);

for n_iter = 2:iteration_max
    rho = NextDensityMatrix(rho, H0matrix, A2_matrix, eigenstate_decomp, dt);
    energy_exp_vec(n_iter) = trace(rho*H0matrix);
    entropy(n_iter) = Entropy(rho);
end

tvector=0:dt:(iteration_max-1)*dt;

figure(2)
plot(tvector,energy_exp_vec,'r-.')
str=['Energy vs t, initial state: ',num2str(bases_a_constrained(ini_ind,:))];
title(str)

figure(3)
plot(tvector,entropy,'b')
str=['Entropy vs t, initial state: ',num2str(bases_a_constrained(ini_ind,:))];
title(str)

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




function energy_deg = Degeneracy_spectra(energies, dim)
%to find the degeneracy of eigen-energies
%input: consecutive differences of eigen-energies, which is a column vector
%output: a two-row matrix with the first row being the energies in the
%spectra and the second row the corresponding degeneracy
global tol_energy;
energy_deg = [energies(1);1];
col = 1;
i=2;
while i<=dim
    if abs(energies(i) - energy_deg(1,col))>tol_energy
        energy_deg = [energy_deg,[energies(i);1]];
        col = col + 1;
    else
        energy_deg(2,col) = energy_deg(2,col) + 1;
    end
    i = i+1;
end
end


function bases_a_constrained = Bases_a_constrained(bases_a, N_sites)
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



function bases_a = Bases_a(N_sites)
%to find the bases in subsystem A
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





function [bases_a_lup, bases_a_ldown] = Symmetry_Zl(bases_a,N_sites)
%to divide the bases into two group based on the spin on the site L
bases_a_lup = [];
bases_a_ldown = [];
total = size(bases_a,1);
for i=1:total
    basis_i = bases_a(i,:);
    if basis_i(N_sites)==1
        bases_a_lup = [bases_a_lup;basis_i];
    else
        bases_a_ldown = [bases_a_ldown;basis_i];
    end
end
end



function H_a_p = Hamiltonian_a_prime(bases,N_sites)
% to find the matrix form of Hamiltonian_a_prime in the bases
%functions used: PXP_prime
dimH = size(bases,1);
H_a_p = zeros(dimH,dimH);
for i=1:dimH
    for j=i+1:dimH
        basis_i = bases(i,:);
        basis_j = bases(j,:);
        H_a_p(j,i)=PXP_prime(basis_i, basis_j, N_sites);
    end
end
H_a_p = H_a_p + H_a_p';    
end

function Hji = PXP_prime(basis_i,basis_j,N_sites)
% the matrix element between basis_i and basis_j
% P_lp1_therm P_{L-1}X_{L} is included in the subsystem A
% used by Hamiltonian_a_prime
global P_lp1_therm;
diff = basis_i - basis_j;
diff_elements = find(diff);

Hji = 0;
if size(diff_elements,2)==1 
    if diff_elements(1)==1
        if basis_i(2)==-1
            Hji = 1;
        end
    elseif diff_elements(1)==N_sites
        if basis_i(N_sites-1)==-1
            Hji = P_lp1_therm;
        end       
    else
        xi = diff_elements(1);
        if basis_i(xi-1)==-1 && basis_i(xi+1)==-1
            Hji = 1;
        end
    end
end
end


function H_a = Hamiltonian_a(bases,N_sites)
% to find the matrix form of Hamiltonian_a in the bases
% functions used: PXP
dimH = size(bases,1);
H_a = zeros(dimH,dimH);
for i=1:dimH
    for j=i+1:dimH
        basis_i = bases(i,:);
        basis_j = bases(j,:);
        H_a(j,i) = PXP(basis_i,basis_j,N_sites);
    end
end
H_a = H_a + H_a';
end


function Hji = PXP(basis_i,basis_j,N_sites)
% the matrix element between basis_i and basis_j
% P_{L-1}X_{L} is excluded in the subsystem A
% used by Hamiltonian_a
diff = basis_i - basis_j;
diff_elements = find(diff);

Hji = 0;
if size(diff_elements,2)==1 
    if diff_elements(1)==1
        if basis_i(2)==-1
            Hji = 1;
        end
    else 
        xi = diff_elements(1);
        if xi < N_sites
            if basis_i(xi-1)==-1 && basis_i(xi+1)==-1
                Hji = 1;
            end
        end
    end
end
end


function A1_omega_matrix = A1_omega_between_eigenstates(eigenstates,energies,bases,N_sites,omega)
% A1_omega = \sum_{a,b} \delta(\omega - (energy_b - energy_a)) |a><a|A1|b><b|
% A1 or A1_omega has no diagonal entries
% functions used: A1_between_eigenstates
tol = 10^(-10);
dim = size(eigenstates,1);
A1_omega_matrix = zeros(dim,dim);
for a=1:dim
    for b = a+1:dim
        if abs(omega - (energies(b) - energies(a)))<tol
            A1_omega_matrix(a,b) = A1_between_eigenstates(eigenstates,bases,b,a,N_sites);
        end
    end
end
A1_omega_matrix = A1_omega_matrix + A1_omega_matrix';
end




function a1_ji = A1_between_eigenstates(eigenstates,bases,i,j,N_sites)
% A_1 = P_{L-1}X_{L}
% to find the matrix element between two 'naturan' bases (\sigma_z) in the real space
% <j|A1|i>
% functions used: A1_between_real_space_bases
% used by: A1_omega_between_eigenstates
tol = 10^(-10);
eigenstate_i = eigenstates(:,i);
eigenstate_j = eigenstates(:,j);
dim = size(bases,1);
a1_ji = 0;
%eigenstates
for n=1:dim %n is the index for the 'natural' bases(\sigma_z) in real space
    eigeni_n = eigenstate_i(n);
    %disp(['n=',num2str(n),' ',num2str(eigeni_n)])
    if abs(eigeni_n) <= tol % eigeni_n==0
        continue
    else
        for m=1:dim %m is the index for the 'natural' bases(\sigma_z) in real space
            eigenj_m = eigenstate_j(m);
            %disp(['m=',num2str(m),' ',num2str(eigenj_m)])
            if m==n || abs(eigenj_m) <= tol %eigenj_m==0
                continue
            else
                a1_ji = a1_ji + eigeni_n*eigenj_m*A1_between_real_space_bases(bases,n,m,N_sites);
            end
        end
    end
end
end



function a1_ji = A1_between_real_space_bases(bases,i,j,N_sites)
% A_1 = P_{L-1}X_{L}
% to find the matrix element between two 'naturan' bases (\sigma_z) in the real space
% <j|A1|i>
% used by A1_between_eigenstates
basis_i = bases(i,:);
basis_j = bases(j,:);
a1_ji = 0;
diff = basis_i - basis_j;
diff_elements = find(diff);
if size(diff_elements,2)==1 
    if diff_elements(1)==N_sites
        if basis_i(N_sites-1)==-1
            a1_ji = 1;
        end
    end
end
end



function eigenstate_decomp = Decompose_eigenpairs(energies, dim)
global tol_energy;
eigenstate_decomp = eigenstate_energy_decomp_positive;
eigenstate_decomp.energies = [];
eigenstate_decomp.eigen_pairs = {};
for i=1:dim
    for j=i+1:dim
        e_diff = energies(j) - energies(i);
        ef = find(abs(eigenstate_decomp.energies - e_diff)<tol_energy);
        if size(ef,2)>0.5 %e_diff in eigenstate_decomp.energies
            ind = ef(1);
            eigenstate_decomp.eigen_pairs{ind} = [eigenstate_decomp.eigen_pairs{ind};[i,j]];
        else
            eigenstate_decomp.energies = [eigenstate_decomp.energies, e_diff];
            eigenstate_decomp.eigen_pairs = [eigenstate_decomp.eigen_pairs, [i,j]];            
        end
    end
end
end




function A2_matrix = A2_Matrix(eigenstates,bases,N_sites)
%to find the A2 matrix in the bases of eigenstates of Hamiltonian
%function uses: A2_between_eigenstates
dim = size(bases,1);
A2_matrix = zeros(dim,dim);
for i=1:dim
    A2_matrix(i,i) = A2_between_eigenstates(eigenstates,bases,i,i,N_sites);
end

for i=1:dim %row index
    for j=i+1:dim %column index
        A2_matrix(i,j) = A2_between_eigenstates(eigenstates,bases,j,i,N_sites);
        A2_matrix(j,i) = conj(A2_matrix(i,j));
    end
end
end



function a2_ji = A2_between_eigenstates(eigenstates,bases,i,j,N_sites)
% A_2 = P_{L}
% to find the matrix element between two 'natural' bases (\sigma_z) in the real space
% <j|A2|i>
% functions used: A2_between_real_space_bases
tol = 10^(-10);
a2_ji = 0;

eigenstate_i = eigenstates(:,i);
eigenstate_j = eigenstates(:,j);
dim = size(bases,1);
%eigenstates
for n=1:dim %n is the index for the 'natural' bases(\sigma_z) in real space
    eigeni_n = eigenstate_i(n);
    %disp(['n=',num2str(n),' ',num2str(eigeni_n)])
    if abs(eigeni_n) <= tol % eigeni_n==0
        continue
    else
        for m=1:dim %m is the index for the 'natural' bases(\sigma_z) in real space
            eigenj_m = eigenstate_j(m);
            %disp(['m=',num2str(m),' ',num2str(eigenj_m)])
            if abs(eigenj_m) <= tol %eigenj_m==0
                continue
            else
                a2_ji = a2_ji + eigeni_n.*conj(eigenj_m).*A2_between_real_space_bases(bases,n,m,N_sites);
            end
        end
    end
end

end



function a2_ji = A2_between_real_space_bases(bases,i,j,N_sites)
% A_2 = P_{L}
% to find the matrix element between two 'naturan' bases (\sigma_z) in the real space
% only diagonal elements survive
% <j|A2|i>
% used by A2_between_eigenstates_brutal_force_method
a2_ji = 0;
if i == j
    if bases(i,N_sites) == -1
        a2_ji = 1;
    end
end
end



function densityMatrix_next = NextDensityMatrix(densityMatrix, H0matrix, A2_matrix, eigenstate_decomp, dt)
%To find density matrix of next step using Runge-Kutta algorithm
%functions used: LindbladEqn
K1=LindbladEqn(densityMatrix, H0matrix, A2_matrix, eigenstate_decomp);
K2=LindbladEqn(densityMatrix+0.5*dt*K1, H0matrix, A2_matrix, eigenstate_decomp);
K3=LindbladEqn(densityMatrix+0.5*dt*K2, H0matrix, A2_matrix, eigenstate_decomp);
K4=LindbladEqn(densityMatrix+dt*K3, H0matrix, A2_matrix, eigenstate_decomp);
densityMatrix_next=densityMatrix + (dt/6)*(K1+2*K2+2*K3+K4);
end



function fmatrix = LindbladEqn(densityMatrix, H0matrix, A2_matrix, eigenstate_decomp)
%The right-hand side(RHS) of the quantum master equation
%only A2
%\hbar=1
%used by NextDensityMatrix
tol = 1e-7;
global imagunit;
dim = size(H0matrix,1);
fmatrix = -imagunit*Commutator(H0matrix,densityMatrix);
dissipator = 0;
n_energy_decomp = size(eigenstate_decomp.energies,2);


for i=1:n_energy_decomp
    e = eigenstate_decomp.energies(i);
    gamma=Gamma(e);
    A2_up = zeros(dim,dim);
    pairs = eigenstate_decomp.eigen_pairs{i};
    ne = size(pairs,1);
    flag=0;
    for j=1:ne
        col_ind = pairs(j,1);
        row_ind = pairs(j,2);
        if A2_matrix(row_ind, col_ind)>tol
            flag=1;
            A2_up(row_ind, col_ind) = A2_matrix(row_ind, col_ind);
        end
    end
    
    if flag>0.5
        A2_down = A2_up';
        dissipator = dissipator + gamma(1)*(A2_up*densityMatrix*(A2_up')-0.5*antiCommutator((A2_up')*A2_up,densityMatrix));
        dissipator = dissipator + gamma(2)*(A2_down*densityMatrix*(A2_down')-0.5*antiCommutator((A2_down')*A2_down,densityMatrix));
    end
end

fmatrix = fmatrix + dissipator;
end



function matrix = Commutator(A,B)
matrix=A*B-B*A;
end

function matrix = antiCommutator(A,B)
matrix=A*B+B*A;
end



function n_b=BoseFun(n)
%n:w=n*energy units
global T
n_b=1./(exp(n./T)-1);
end


function Gamma_n=Gamma(n)
global gamma
n_b=BoseFun(n);
Gamma1=gamma.*n_b; %system go up in energy
Gamma2=gamma.*(n_b+1); %system go down in energy
Gamma_n=[Gamma1,Gamma2];
end


