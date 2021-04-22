function level_statistics1
%periodic boundary condition
%fermions on two chains
%Nparticles1 <= Nparticles2
%With inversion symetry: T1_hop=T2_hop
%account for inversion symmetry

clear all

%%parameters
global Nsites Nparticles1 Nparticles2 T1_hop T2_hop
Nsites=6;
Nparticles1=3;
Nparticles2=3;
T1_hop=-1;
T2_hop=T1_hop;


%%sectors
sectors=[];
if Nparticles1==0||Nparticles1==1
    sectors=[sectors,Nparticles2];
else
    seed_sec=integer_partitions(Nparticles2,Nparticles1)
    for i=1:size(seed_sec,1)
        seedi=seed_sec(i,:);
        seedi_perm = perms(seedi(2:end));
        seedi_perm = [ repmat(seedi(1),size(seedi_perm,1),1) seedi_perm ];
        sectorsi=DeleteDuplicates(seedi_perm);
        sectors=[sectors;sectorsi];
    end
end
sectors 



%%state matrix for each sector



fignum=1;
for num_sec=1:size(sectors,1)
pause
sector=sectors(num_sec,:)
bases_sector=States_sector(1, sector);
trans_rule=InverseTrans_bases(bases_sector);
[self_even,pair_inv]=Splitbases_Inversion(trans_rule)

H_spin1=Hamiltonian_spin1(bases_sector);
H_spin2=Hamiltonian_spin2(bases_sector);
H=H_spin1+H_spin2;

%%how many nonzero elements in each column of H
Num_nonzero_H=zeros(1,size(H,1));
for i=1:size(H,1)
    col_Hi=H(:,i);
    Num_nonzero_H(i)=size(find(col_Hi==T1_hop|col_Hi==T2_hop),1);
end
%Num_nonzero_H
%size(nonzeros(diff(Num_nonzero_H)),1)==0

Num_nonzero_H_spin1=zeros(1,size(H_spin1,1));
for i=1:size(H_spin1,1)
    col_Hi=H_spin1(:,i);
    Num_nonzero_H_spin1(i)=size(find(col_Hi==T1_hop),1);
end
%Num_nonzero_H_spin1
%size(nonzeros(diff(Num_nonzero_H_spin1)),1)==0

Num_nonzero_H_spin2=zeros(1,size(H_spin2,1));
for i=1:size(H_spin2,1)
    col_Hi=H_spin2(:,i);
    Num_nonzero_H_spin2(i)=size(find(col_Hi==T2_hop),1);
end
%Num_nonzero_H_spin2
%size(nonzeros(diff(Num_nonzero_H_spin2)),1)==0



%%nonzeros_H: column i tells what states the ith state in the bases can be
%%mapped to under H
nonzeros_H=zeros(max(Num_nonzero_H),size(H,1));
for i=1:size(H,1)
    col_Hi=H(:,i);
    nonzero_H_i=find(col_Hi==T1_hop);
    nonzeros_H(:,i)=[nonzero_H_i;zeros(max(Num_nonzero_H)-size(nonzero_H_i,1),1)];
end
%nonzeros_H


nonzeros_H_spin1=zeros(max(Num_nonzero_H_spin1),size(H_spin1,1));
for i=1:size(H_spin1,1)
    col_Hi=H_spin1(:,i);
    nonzero_H1_i=find(col_Hi==T1_hop);
    nonzeros_H_spin1(:,i)=[nonzero_H1_i;zeros(max(Num_nonzero_H_spin1)-size(nonzero_H1_i,1),1)];
end
%nonzeros_H_spin1

nonzeros_H_spin2=zeros(max(Num_nonzero_H_spin2),size(H_spin2,1));
for i=1:size(H_spin2,1)
    col_Hi=H_spin2(:,i);
    nonzero_H2_i=find(col_Hi==T2_hop);
    nonzeros_H_spin2(:,i)=[nonzero_H2_i;zeros(max(Num_nonzero_H_spin2)-size(nonzero_H2_i,1),1)];
end
%nonzeros_H_spin2


H_even=Hevensector(self_even,pair_inv,nonzeros_H);
H_odd=Hoddsector(pair_inv,nonzeros_H);


% [eigenstates,energies_diag]=eig(H); %eigenstates are column vectors


[eigenstates_even,energy_even]=eig(H_even);
% format long 
%transpose(diag(energy_even))
%eigenstates_even
e_even_diff=diff(diag(energy_even));

[eigenstates_odd,energy_odd]=eig(H_odd);
% format long
%transpose(diag(energy_odd))
%eigenstates_odd
e_odd_diff=diff(diag(energy_odd));


figure(fignum)
fignum=fignum+1;
hi_even=histogram(diag(energy_even));
hi_even.BinWidth=min(nonzeros(diag(e_even_diff)));
title(['degeneracy for the even sector of ',num2str(sector), ' (', num2str([Nsites, Nparticles1, Nparticles2]),')'])
%title(['level statistics for the even sector ',num2str(sector), ' (', num2str([Nsites, Nparticles1, Nparticles2]),')'])

figure(fignum)
fignum=fignum+1;
hi_odd=histogram(diag(energy_odd));
hi_odd.BinWidth=min(nonzeros(diag(e_odd_diff)));
title(['degeneracy for the odd sector of ',num2str(sector), ' (', num2str([Nsites, Nparticles1, Nparticles2]),')'])
%title(['level statistics for the odd sector ',num2str(sector), ' (', num2str([Nsites, Nparticles1, Nparticles2]),')'])

end



% for i=1:size(eigenstates,2)
%     energies_diag(i,i)
%     OccupationNumber(bases_sector, eigenstates(:,i))
% end



end


function H_even=Hevensector(self_even,pair_inv,nonzeros_H)
global T1_hop
size_self=size(self_even,2);
size_pair=size(pair_inv,2);
dimH=size_self+size_pair;
H_even=zeros(dimH,dimH);
selfUpair= [[self_even;zeros(1,size(self_even,2))],pair_inv];
for i=1:dimH
    basis_i=selfUpair(:,i);
    rep_i=basis_i(1);
    col_rep=nonzeros_H(:,rep_i);
    for k=col_rep'
        if k==0
            break
        end
        [row,j]=find(selfUpair==k);
        if (i<=size_self && j<=size_self) || (i>size_self && j>size_self)
            H_even(i,j)=T1_hop;
        else
            H_even(i,j)=sqrt(2)*T1_hop;
        end
    end
end
if H_even~=transpose(H_even)
    disp('H_even is a Hermitian matrix!')
end
end


function H_odd=Hoddsector(pair_inv,nonzeros_H)
global T1_hop
dimH=size(pair_inv,2);
H_odd=zeros(dimH,dimH);
for i=1:dimH
    basis_i=pair_inv(:,i);
    rep_i=basis_i(1);
    col_rep=nonzeros_H(:,rep_i);
    for k=col_rep'
        if k==0
            break
        end
        [row,j]=find(pair_inv==k);
        H_odd(i,j)=T1_hop;
    end
end
if H_odd~=transpose(H_odd)
    disp('H_odd is a Hermitian matrix!')
end
end


function [self_even,pair_inv]=Splitbases_Inversion(trans_rule)
self_even=[];
pair_inv=[];
for i=1:size(trans_rule,2)
    if trans_rule(1,i)==trans_rule(2,i)
        self_even=[self_even,trans_rule(1,i)];
    else
        pair_inv=[pair_inv,trans_rule(:,i)];
    end
end
pair_inv=deleteDupPairs(pair_inv);
end

function pairs2=deleteDupPairs(pairs1)
pairs2=[];
pairs2=[pairs2,pairs1(:,1)];
for i=2:size(pairs1,2)
    pairs1_i=pairs1(:,i);
    flag=1;
    for j=1:size(pairs2,2)
        if sum(pairs2(:,j)==flip(pairs1_i))==2
            flag=0;
            break
        end
    end
    if flag==1
        pairs2=[pairs2,pairs1_i];
    end
end

end


function occupationNum=OccupationNumber(bases, eigenstate)
%occupation number for each site for the given eigen state
global Nsites
occupationNum=zeros(2,Nsites);
A2=eigenstate.*eigenstate; %col vector
for i=1:Nsites
    ni_spin2=bases(1:2:size(bases,1)-1,i);%col vector
    ni_spin1=bases(2:2:size(bases,1),i); %col vector
    occupationNum(1,i)=A2'*ni_spin2;
    occupationNum(2,i)=A2'*ni_spin1;
end
end


function flag=IsEqualEigenStates(v1,v2)
%eigen vectors are column vectors
tol=0.0001;
diff=abs(v1-v2);
if max(diff)<tol
    flag=1;
else
    flag=0;
end
end

function eigenstate_trans=InverseTrans_eigenstates(eigenstate,trans_rule)
%Inversion transformation of eigenstates on a certain bases
n_eigen=size(eigenstate,1);
eigenstate_trans=zeros(n_eigen,1);
for i=1:n_eigen
    eigenstate_trans(trans_rule(2,i))=eigenstate(i);
end
end



function trans_rule=InverseTrans_bases(bases)
global Nsites
%Inversion trasnformation of bases
nstates=size(bases,1)/2;
trans_rule=zeros(2,nstates);
trans_rule(1,:)=1:nstates;
trans_rule2=zeros(1,nstates);
for i=1:nstates
    if ismember(i,trans_rule2)
        trans_rule2(i)=find(trans_rule2==i);
        continue
    end
    
    basis_i=bases(2*(i-1)+1:2*i,:);
    basis_i_trans=InverseTrans(basis_i);
    for j=1:nstates
        basis_j=bases(2*(j-1)+1:2*j,:);
        if sum(sum(basis_i_trans==basis_j))==2*Nsites
            trans_rule2(i)=j;
            break
        end
    end
end
trans_rule(2,:)=trans_rule2;
end


function state_f=InverseTrans(state_i)
%Inversion transformation of one state 
%site j --> N-j+1
%species (1) <--> (2)
state_f_spin2=state_i(2,:);
state_f_spin1=state_i(1,:);
state_f_spin1=flip(state_f_spin1);
state_f_spin2=flip(state_f_spin2);
state_f=[state_f_spin2;state_f_spin1];
end


function H=Hamiltonian_spin1(stateMatrix)
%<row| H | col>

dimH=size(stateMatrix,1)/2;
H1=zeros(dimH,dimH);
for row=1:dimH
    state_r_spin2=stateMatrix(2*(row-1)+1,:);
    state_r_spin1=stateMatrix(2*row,:);
    for col=row+1:dimH
%         state_r_spin2=stateMatrix(2*(row-1)+1,:)
%         state_r_spin1=stateMatrix(2*row,:)
        state_c_spin2=stateMatrix(2*(col-1)+1,:);
        state_c_spin1=stateMatrix(2*col,:);
%         pause
        H1(row,col)=H_spin1Hop(state_r_spin1,state_r_spin2,state_c_spin1,state_c_spin2);
    end
end
H2=transpose(H1);
H=H1+H2;
end


function H=Hamiltonian_spin2(stateMatrix)
%<row| H | col>

dimH=size(stateMatrix,1)/2;
H1=zeros(dimH,dimH);
for row=1:dimH
    state_r_spin2=stateMatrix(2*(row-1)+1,:);
    state_r_spin1=stateMatrix(2*row,:);
    for col=row+1:dimH
%         state_r_spin2=stateMatrix(2*(row-1)+1,:)
%         state_r_spin1=stateMatrix(2*row,:)
        state_c_spin2=stateMatrix(2*(col-1)+1,:);
        state_c_spin1=stateMatrix(2*col,:);
%         pause
        H1(row,col)= H_spin2Hop(state_r_spin1,state_r_spin2,state_c_spin1,state_c_spin2);
    end
end
H2=transpose(H1);
H=H1+H2;
end



function Hspin1hop_rc=H_spin1Hop(state_r_spin1,state_r_spin2,state_c_spin1,state_c_spin2)
global T1_hop Nsites
Hspin1hop_rc=0;
diff_spin2=state_r_spin2-state_c_spin2;
diff_spin1=state_r_spin1-state_c_spin1;
pos_p=find(diff_spin1==1);
pos_h=find(diff_spin1==-1);
% sum(abs(diff_spin2))==0
% size(pos_p,2)==1
% size(pos_h,2)==1
% mod(min([pos_h-1,pos_p-1])+1,Nsites)==max([pos_h-1,pos_p-1])
if sum(abs(diff_spin2))==0 && size(pos_p,2)==1 && size(pos_h,2)==1 
    if (mod(pos_h,Nsites)==(pos_p-1)  && state_c_spin2(pos_h)==0) || (mod(pos_p,Nsites)==(pos_h-1) && state_c_spin2(pos_p)==0)
        Hspin1hop_rc=T1_hop;
%         pause
    end
end
end

function Hspin2hop_rc=H_spin2Hop(state_r_spin1,state_r_spin2,state_c_spin1,state_c_spin2)
global T2_hop Nsites
Hspin2hop_rc=0;
diff_spin2=state_r_spin2-state_c_spin2;
diff_spin1=state_r_spin1-state_c_spin1;
pos_p=find(diff_spin2==1);
pos_h=find(diff_spin2==-1);
% sum(abs(diff_spin1))==0
% size(pos_p,2)==1
% size(pos_h,2)==1
% mod(min([pos_h-1,pos_p-1])+1,Nsites)==max([pos_h-1,pos_p-1])
if sum(abs(diff_spin1))==0 && size(pos_p,2)==1 && size(pos_h,2)==1 
    if (mod(pos_h,Nsites)==(pos_p-1) && state_c_spin1(pos_p)==0) || (mod(pos_p,Nsites)==(pos_h-1) && state_c_spin1(pos_h)==0)
        Hspin2hop_rc=T2_hop;
%         pause
    end
end
end


function sectorsi=DeleteDuplicates(seedi_perm)
%select sectors up to cyclic permutation
if size(seedi_perm,1)==1
    sectorsi=seedi_perm;
else
    sectorsi=[];
    seedi_perm=flip(seedi_perm);
    original=seedi_perm(1,:);
    sectorsi=[sectorsi;original];
    perms=seedi_perm(2:end,:);
    for i=1:size(perms,1)
        if ~IsCyclicEqual(sectorsi,perms(i,:))
            sectorsi=[sectorsi;perms(i,:)];
        end
    end
end
end


function flag=IsCyclicEqual(states1, state2)
if size(states1,2)~=size(state2,2)
    disp('error');
else
    flag=0;
    for j=1:size(states1,1)
        states1_j=states1(j,:);
        for i=0:1:(size(states1,2)-1)
            if sum(states1_j==circshift(state2,i))==size(states1,2)
                flag=1;
                break
            end
        end
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








function stateMatrix=States_sector(indspin1_start, sector)
% sum(sector)=Nparticles2
% length(sector)=Nparticles1
% every two rows represent one state (bot row--spin 1;top row--spin 2)
global Nsites Nparticles1 

primeNums=[5,7,11];
%for a given sector, how choices for the indices of spin 1
len=Nparticles1;
sector_replaced=sector;
sector_replaced(sector_replaced==0)=1;
Xrange=cell(1,len);
for i=1:len
    sectorR_new=sector_replaced;
    sectorR_new(i)=[];
    inti=sector_replaced(i);
    endi=Nsites-sum(sectorR_new);
    rangei=inti:endi;
    Xrange{i}=rangei;
end
[A{1:len}]=ndgrid(Xrange{:});
B = reshape(cat(len+1,A{:}),[],len);
indC = sum(B,2)==Nsites;

% disp('index of spin1')
indspin1_seed=B(indC,:);
% pause

indspin1 = DeleteIndspin1Duplicates(sector,indspin1_seed);
% pause

size_indspin1=size(indspin1);
nstates_spin1=size_indspin1(1);

maxSteps_circshift=zeros(1,nstates_spin1);
for i=1:nstates_spin1
    maxSteps_circshift(i)=MaxSteps_circshift(sector,indspin1(i,:));
end


states_spin1=zeros(nstates_spin1,Nsites);
for i=1:nstates_spin1
    statei_spin1=zeros(1,Nsites);
    indspin1_i=indspin1(i,:);
    statei_spin1(indspin1_start)=1;
    indspin1_j=indspin1_start;
    for j=1:len-1
        indspin1_j=AddInd(indspin1_j,indspin1_i(j));
        statei_spin1(indspin1_j)=1;               
    end
    states_spin1(i,:)=statei_spin1;
end

% disp('states of spin1')
% states_spin1
% pause



%choices of indices of spin 2
stateMatrix=[];
for i=1:nstates_spin1
%     disp(['state ',num2str(i),' of spin1'])
%     states_spin1(i,:);
%     disp(['index of state ',num2str(i),' of spin1'])
    indspin1_i=indspin1(i,:);
%     pause
    
    indspin2_choices=cell(1,len);
    indspin1_startj=indspin1_start;
    for j=1:len
        choicesj=0:1:indspin1_i(j)-1;
        indspin2_choices{j}=AddInd(choicesj,indspin1_startj);
        indspin1_startj=indspin1_startj+indspin1_i(j);
    end
    
    indspin2={};
    for j=1:len
        nspin2_j=sector(j);
        if nspin2_j~=0
            indspin2=[indspin2, nchoosek(indspin2_choices{j},nspin2_j)];
        else
            continue              
        end
    end
    
%     disp('index of spin2')
%     for j=1:size(indspin2,2)
%         indspin2{j}
%     end
%     pause
    
    col_indspin2=size(indspin2,2);
    rind=cell(1,col_indspin2);
    for j=1:col_indspin2
        rindj=size(indspin2{j},1);
        rind{j}=1:rindj;
    end
    [R{1:col_indspin2}]= ndgrid(rind{:});
    R1=reshape(cat(col_indspin2+1,R{:}),[],col_indspin2);
    
    indspin2_cated=[];
    for j=1:size(R1,1)
        indspin2_catedj=[];
        for k=1:col_indspin2
            indspin2_catedj=[indspin2_catedj, indspin2{k}(R1(j,k) ,:)];
        end
        indspin2_cated=[indspin2_cated;indspin2_catedj];
    end
%     disp('index of spin2 (concatenated)')
%     indspin2_cated
    
    states_spin2_i=zeros(size(indspin2_cated,1),Nsites);
    for j=1:size(indspin2_cated,1)
        for k=1:size(indspin2_cated,2)
            states_spin2_i(j, indspin2_cated(j,k))=1;
        end
    end
%     disp('states of spin2')
%     states_spin2_i
    
    
    states_combined_1=CombineStates(states_spin1(i,:), states_spin2_i);
    states_combined_allcircshifted=[];
    states_combined_allcircshifted=[states_combined_allcircshifted;states_combined_1];
    for shiftstep=1:maxSteps_circshift(i)-1
        states_combined_allcircshifted=[states_combined_allcircshifted;  circshift(states_combined_1,shiftstep,2)];
    end
    
    stateMatrix=[stateMatrix;states_combined_allcircshifted];
end
end



function maxStep_circshift=MaxSteps_circshift(sector,indspin1)
indspin1_matrix=[sector; indspin1];
for i=1:size(indspin1,2)
    if sum(sum(indspin1_matrix==circshift(indspin1_matrix,i,2)))==2*size(indspin1,2)
        maxStep_circshift=sum(indspin1(1:i));
        break
    end
end
end



function indspin1 = DeleteIndspin1Duplicates(sector,indspin1_seed)

n_seed=size(indspin1_seed,1);
indspin1=[];
indspin1=[indspin1;indspin1_seed(1,:)];
for i=2:n_seed
    indspin1_seedi=indspin1_seed(i,:);
    if ~IsEqualIndspin1(indspin1,indspin1_seedi,sector)
        indspin1=[indspin1;indspin1_seedi];
    end
end
end


function flag=IsEqualIndspin1(indspin1,indspin1_seedi,sector)
flag=0;
size_ind=size(indspin1);
row=size_ind(1);
col=size_ind(2);
indspin1_seedi_matrix=[sector;indspin1_seedi];
for i=1:row
    indspin1_i=indspin1(i,:);
    indspin1_i_matrix=[sector;indspin1_i];
    for j=0:col-1
        if sum(sum(indspin1_seedi_matrix==circshift(indspin1_i_matrix,j,2)))==2*col
            flag=1;
            break
        end
    end
    if flag==1
        break
    end
end

end


function statesCombined=CombineStates(state_spin1, states_spin2)
global Nsites
nr=size(states_spin2,1);
statesCombined=zeros(2*nr,Nsites);
for i=1:nr
    statesCombined((i-1)*2+1:i*2,:)=[states_spin2(i,:);state_spin1];
end

end


function ind=AddInd(ind1, ind2)
global Nsites
ind=mod(ind1+ind2,Nsites);
ind(ind==0)=Nsites;
end
