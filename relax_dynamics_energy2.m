function relax_dynamics_energy2
% Relaxation dynamics with only 1 channel for g=1/2, g=2/3 and spinless
% fermion

clear all

global T gamma imagunit dt dimH
imagunit=sqrt(-1);
T=0.01; %temperature
gamma=10;%gamma=2\pi*Density of states*(coupling constant)^2
dt=0.001;


%%parameters

Eint=5; %Eint>0;  initial energy with the convention that the energy of the groud state is 0. Eint>0
nmax=1000;
channels=[2]; % the channel openned, which is  the number of energy units partilce can go up or down in energy space
size_channels=size(channels);
nchannel=size_channels(2);


%%state matrix for Haldane statistics particles, which can be represented by integer
%%partitions of energies pumped into the systems
[states_h,deg]=partitionMatrix(Eint);
states_h
deg
%k_int = 17;
%states_h(k_ini,:)

%%state matrix for spinless fermions
size_belowfl=Eint+1;
size_abovefl=Eint;
groundState_f=[ones(1,size_belowfl),zeros(1,size_abovefl)];
states_f=stateMatrix_fermi(groundState_f,states_h,size_belowfl)



%%Hamiltonian of the system, H0
dimH=sum(deg);
H0=zeros(dimH, dimH);
size_deg=size(deg);
size_degc=size_deg(2);
for j=2:size_degc
    count_start=sum(deg(1:j-1));
    for k=1:deg(j)
        H0(count_start+k,count_start+k)=j-1;
    end
end
H0


%Gammaton
Gamma_matrix=zeros(nchannel,2);
for i=1:nchannel
    channel=channels(i);
    Gamma_matrix(i,:)=Gamma(channel);
end
disp("gamma matrix:")
disp(Gamma_matrix)


%%Lindblad operator 1 for Haldane statistics particles
L1ac_h=zeros(dimH,dimH,nchannel);
for i=1:nchannel
    channel=channels(i);
    L1ac_h(:,:,i)=Lindblad_h(states_h,deg,channel);
end
disp("L1 of HES:")
L1ac_h



%%Lindblad operator 1 for spinless fermions
L1ac_f=zeros(dimH,dimH,nchannel);
for i=1:nchannel
    channel=channels(i);
    L1ac_f(:,:,i)=Lindblad_f(states_f,deg,channel);
end
disp("L1 of spinless fermion:")
L1ac_f


%%initial density matrix 
Rho_h=zeros(dimH,dimH);
Rho_f=zeros(dimH,dimH);
k_ini = 16;
states_h(k_ini,:)
Rho_h(k_ini,k_ini)=1;
Rho_f(k_ini,k_ini)=1;
% subspace =[17,18,19];
% dim_subspace = size(subspace, 2);
%generate random numbers
% diag_rand = zeros(1,dim_subspace);
% for kk_ind=1:dim_subspace
%     diag_rand(kk_ind) = random('Uniform',0,1);
% end
% total_diag_rand = sum(diag_rand);
% for kk_ind=1:dim_subspace
%     diag_rand(kk_ind) = diag_rand(kk_ind)./total_diag_rand;
% end
% offdiag_rand_re = [];
% offdiag_rand_imag = [];
% for kr_ind=1:dim_subspace
%     for kc_ind=kr_ind+1:dim_subspace
%         offdiag_rand_re = [offdiag_rand_re, random('Uniform',0,1)];
%         offdiag_rand_imag = [offdiag_rand_imag, random('Uniform',0,1)];
%     end
% end
% disp("random diag:")
% fprintf('%10.9f ', diag_rand)
% fprintf('\n')
% disp("random offdiag re:")
% fprintf('%10.9f ', offdiag_rand_re)
% fprintf('\n')
% disp("random offdiag imag:")
% fprintf('%10.9f ', offdiag_rand_imag)
% fprintf('\n')

% diag_rand = [0.329187283 0.260412561 0.410400156];
% offdiag_rand_re = [0.627973359 0.431651170 0.984063724];
% offdiag_rand_imag = [0.291984080 0.015487126 0.167168410];
% for kk_ind=1:dim_subspace
%     kk = subspace(kk_ind);
%     Rho_f(kk,kk)=diag_rand(kk_ind);
%     Rho_h(kk,kk)=Rho_f(kk,kk);
% end
% count_offdiag = 0;
% for kr_ind=1 : dim_subspace
%     kr = subspace(kr_ind);
%     for kc_ind=kr_ind+1 : dim_subspace
%         kc = subspace(kc_ind);
%         count_offdiag = count_offdiag + 1;
%         re = offdiag_rand_re(count_offdiag);
%         imag = offdiag_rand_imag(count_offdiag);
%         Rho_f(kr,kc) = re + imag.*imagunit;
%         Rho_h(kr,kc) = Rho_f(kr,kc);
%         Rho_f(kc,kr) = re - imag.*imagunit;
%         Rho_h(kc,kr) = Rho_f(kc,kr);
%     end
% end
% disp("initial density matrix of spinless fermions/semions:")
% Rho_f(subspace,subspace)
% Rho_h(subspace,subspace)
% for i=1:10
%     Rho_h(i,i)=0.1;
% end
% for i=1:10
%     Rho_f(i,i)=0.1;
% end



%%diagonal elements of Rho
diag_h=zeros(dimH,nmax);
diag_f=zeros(dimH,nmax);
for j=1:dimH
    diag_h(j,1)=Rho_h(j,j);
    diag_f(j,1)=Rho_f(j,j);
end


%%Energy of  initial state of spinless fermion
energyVector_h=zeros(1,nmax);
energyVector_h(1)=trace(Rho_h*H0);
energyVector_f=zeros(1,nmax);
energyVector_f(1)=trace(Rho_f*H0);

%%Lindbald dynamics (Runge-kutta method)
for i=2:nmax
   Rho_f=NextDensityMatrix(H0,L1ac_f,Gamma_matrix,Rho_f,nchannel);
   Rho_h=NextDensityMatrix(H0,L1ac_h,Gamma_matrix,Rho_h,nchannel);
   
   energyVector_f(i)=trace(Rho_f*H0);
   energyVector_h(i)=trace(Rho_h*H0);
   
   for j=1:dimH
    diag_h(j,i)=Rho_h(j,j);
    diag_f(j,i)=Rho_f(j,j);
   end
end






tvector=0:dt:(nmax-1)*dt;

% kf=1;
%spinless fermion
% V1 = VideoWriter('Haldane statistical particles.avi');
% % V1 = VideoWriter('spinless fermion.avi');
% V1.FrameRate = 20;
% open(V1);
% for i=1:nmax
% %     figure
%     diagi=diag_h(:,i);
%     bar(diagi');
%     xticklabels({'0','1','2','2','3','3','3','4','4','4','4','4'});
%     ylim([0,1]);
%     frame = getframe(gcf);
%     writeVideo(V1,frame);
% end
% close(V1);




% 
% kf=kf+1;
% figure(kf)
% count_start=0;
% for i=1:size_degc
%     degi=deg(i);
%     for j=1:degi
%         k=count_start+j;
%         subplot(3,4,k);
%         plot(tvector,diag_f(k,:),'r-.')
%         axis([0 0.5 0 1])
%         title(['E=',num2str(i-1),',Rho(',num2str(k),',',num2str(k), '), '])
%         hold off
%     end
%     count_start=count_start+degi;
% end





%semion (haldane)
%kf=kf+1;
% kf=1;
% figure(kf)
% count_start=0;
% for i=1:size_degc
%     degi=deg(i);
%     for j=1:degi
%         k=count_start+j;
%         subplot(3,4,k);
%         plot(tvector,diag_h(k,:),'b')
%         axis([0 1.5 0 1])
%         title(['E=',num2str(i-1),',Rho(',num2str(k),',',num2str(k), '), '])
%         hold off
%     end
%     count_start=count_start+degi;
% end




% kf=kf+1;
kf=1;
figure(kf)
plot(tvector,energyVector_f,'r-.')
hold on
plot(tvector,energyVector_h,'b-')
hold on
plot(tvector,energyVector_h,'g--')
title(['energy vs time', ', k=',num2str(channels)])
hold off
legend('Spinless fermion','g=1/2','g=2/3')



end




function densityMatrix_next=NextDensityMatrix(H0matrix,L1ac,Gamma_matrix,densityMatrix,nchannel)
%To find density matrix of next step in Runge-Kutta algorithm
global dt
K1=LindbladEqn(H0matrix,L1ac,Gamma_matrix,densityMatrix,nchannel);
K2=LindbladEqn(H0matrix,L1ac,Gamma_matrix,densityMatrix+0.5*dt*K1,nchannel);
K3=LindbladEqn(H0matrix,L1ac,Gamma_matrix,densityMatrix+0.5*dt*K2,nchannel);
K4=LindbladEqn(H0matrix,L1ac,Gamma_matrix,densityMatrix+dt*K3,nchannel);
densityMatrix_next=densityMatrix + (dt/6)*(K1+2*K2+2*K3+K4);
end


function fmatrix=LindbladEqn(H0matrix,L1ac,Gamma_matrix,densityMatrix,nchannel)
%make the right-hand side of Lindblad equation
global imagunit dimH

fmatrix=-imagunit*Commutator(H0matrix,densityMatrix);
for i=1:nchannel
    L1_i=L1ac(:,:,i);
    Gamma_i=Gamma_matrix(i,:);
    for j=1:dimH-1
        for k=j+1:dimH
            if L1_i(k,j)~=0
                L1_is = zeros(dimH,dimH);
                L1_is(k,j) = L1_i(k,j);
                L2_is = transpose(L1_is);
                %pause
                f_dissipation1 = Gamma_i(1)*(L1_is*densityMatrix*transpose(L1_is)-0.5*antiCommutator(transpose(L1_is)*L1_is,densityMatrix));
                f_dissipation2 = Gamma_i(2)*(L2_is*densityMatrix*transpose(L2_is)-0.5*antiCommutator(transpose(L2_is)*L2_is,densityMatrix));
                fmatrix = fmatrix + f_dissipation1 + f_dissipation2;
                %pause
            end
        end
    end   
end
end



function matrix=Commutator(A,B)
matrix=A*B-B*A;
end

function matrix=antiCommutator(A,B)
matrix=A*B+B*A;
end





function L1_h=Lindblad_h(states_h,deg,kchannel)
%Lindblad operator 1 for haldane statistics particles
dimL=sum(deg);
L1_h=zeros(dimL,dimL);
size_deg=size(deg);
len_deg=size_deg(2);
cstart=0; % low
rstart=sum(deg(1:kchannel)); %high
for j=1:(len_deg-kchannel)
    for k1=1:deg(j)
        for k2=1:deg(j+kchannel)
            col=cstart+k1;
            row=rstart+k2;
            statelow=states_h(col,:);
            statehigh=states_h(row,:);
            L1_h(row,col)=L1TwoStates_h(statelow, statehigh, kchannel);
        end
    end
    cstart=cstart+deg(j);
    rstart=rstart+deg(j+kchannel);
end
end


function element=L1TwoStates_h(statelow, statehigh, kchannel)
diffstate=statehigh-statelow;
nonzero_diff=nonzeros(diffstate);
size_nonzero=size(nonzero_diff);
if size_nonzero(1)==1 && sum(nonzero_diff)==kchannel
    element=1;
else
    element=0;
end
end



function L1_f=Lindblad_f(states_f,deg,kchannel)
%Lindblad operator 1 for spinless fermions
dimL=sum(deg);
L1_f=zeros(dimL,dimL);
size_deg=size(deg);
len_deg=size_deg(2);
cstart=0; % low
rstart=sum(deg(1:kchannel)); %high
for j=1:(len_deg-kchannel)
    for k1=1:deg(j)
        for k2=1:deg(j+kchannel)
            col=cstart+k1;
            row=rstart+k2;
            statelow=states_f(col,:);
            statehigh=states_f(row,:);
            L1_f(row,col)=L1TwoStates_f(statelow, statehigh, kchannel);
        end
    end
    cstart=cstart+deg(j);
    rstart=rstart+deg(j+kchannel);
end
end

function element=L1TwoStates_f(statelow, statehigh, kchannel)
diffstate=statehigh-statelow;
ind_low=find(diffstate==-1);
ind_high=find(diffstate==1);
if (ind_high-ind_low)==kchannel
    n_hopover=sum(statelow(ind_low:ind_high))-1;
    element=(-1)^n_hopover;
else
    element=0;
end
end







function n_b=BoseFun(kchannel)
%n:w=n*energy units
global T
n_b=1./(exp(kchannel./T)-1);
end



function Gamma_n=Gamma(n)
global gamma
n_b=BoseFun(n);
Gamma1=gamma.*n_b;
Gamma2=gamma.*(n_b+1);
Gamma_n=[Gamma1,Gamma2];
end



function states_f=stateMatrix_fermi(groundState_f,stateMatrix_Haldane,size_belowfl)
%output: state matrix for spinless fermions
sizep=size(stateMatrix_Haldane);
states_f=groundState_f;
for j=2:sizep(1)
    Evector=stateMatrix_Haldane(j,:);
    statej=MakeExcitedStates(groundState_f,Evector,size_belowfl);
    states_f=[states_f;statej];
end
end


function state=MakeExcitedStates(groundState_f,Evector,fl)
% output: excited state with an energy assignment(Evector) for each fermion 
% input: fl--fermi level
state=groundState_f;
level=fl;
for e = Evector
    if e==0
        break
    end
    state(level)=0;
    state(level+e)=1;
    level=level-1;
end
end


function [integerParMatrix,deg]=partitionMatrix(n)
%output: the matrix of integer partitions of all number less or equal to n
%input: an integer n
if n==0
    integerParMatrix=integer_partitions(n);
    deg=[1];
else
    Matrix_n=integer_partitions(n);
    size_n=size(Matrix_n);
    [Matrix_nless,deg1]=partitionMatrix(n-1);
    size_nless=size(Matrix_nless);
    deg=[deg1,size_n(1)];
    diff_col=abs(size_n(2)-size_nless(2));
    if diff_col==0
        integerParMatrix=[Matrix_nless;Matrix_n];
    else
        row_nless=size_nless(1);
        zerosAppend=zeros(row_nless,1);
        Matrix_nless_fixed=[Matrix_nless,zerosAppend];
        integerParMatrix=[Matrix_nless_fixed;Matrix_n];
    end
end
end






function S = integer_partitions(n,count)
%integer partitions of n
if nargin == 1
    count = n;
end
if n < 0 || n ~= round(n)
    error('Only nonnegative integers allowed!');
elseif n == 0
    if count == 0
        S = 0;
    else
        S = zeros(1,count);
    end
else
    x = ones(1,n);
    x(1) = n;
    m = 1;
    h = 1;
    M = [x(1:m) zeros(1,n-m)];
    while x(1) ~= 1
        if x(h) == 2 
           m = m + 1;
           x(h) = 1;
           h = h - 1;
        else
           r = x(h) - 1;
           t = m - h + 1;
           x(h) = r;
           while t >= r
               h = h + 1;
               x(h) = r;
               t = t - r;
           end
           if t == 0
               m = h;
           else
               m = h + 1;
               if t > 1
                   h = h + 1;
                   x(h) = t;
               end
           end
        end
        M = cat(1,M,[x(1:m) zeros(1,n-m)]);
    end
    if count > n
        M = cat(2,M,zeros(size(M,1),count-n));
    end
    S = [];
    for i = 1:size(M,1)
        if(sum(M(i,1:count)) == n)
            S = cat(1,S,M(i,1:count));
        end
    end
end
end