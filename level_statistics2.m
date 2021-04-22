function level_statistics2
%periodic boundary condition
%fermions on two chains
%Nparticles1 <= Nparticles2
%with or without inversion symmetry (T1_hop = or~= T2_hop or N1 = or~= N2)
%account for translational symmetry (Total momentum)

clear all

%%parameters
global Nsites Nparticles1 Nparticles2 T1_hop T2_hop imagUnit omega
imagUnit=sqrt(-1);
Nsites=8;
Nparticles1=4;
Nparticles2=4;
T1_hop=-1;
T2_hop=-0.7;
omega=exp(imagUnit*2*pi/Nsites);

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
% sectors 



%%state matrix for each sector
fignum=1;
%sector_considered=sectors;
%sector_considered=seed_sec;
sector_considered=[1,1,1,1];
for i=1:size(sector_considered,1)
    sector=sector_considered(i,:)
    [bases_sector,bases_fourier_sector]=States_sector(1, sector);
    %size(bases_sector,1)/2
    %bases_fourier_sector.num
    %celldisp(bases_fourier_sector.states)
    bases_fourier_sector.momenta=Momenta(bases_fourier_sector,Nsites);
    %celldisp(bases_fourier_sector.momenta)
    
    disp("bases found")
    
    H_spin1=Hamiltonian_spin1(bases_sector);
    H_spin2=Hamiltonian_spin2(bases_sector);
%     H=H_spin1+H_spin2;
%     Num_nonzero_H=zeros(1,size(H,1));
%     for j=1:size(H,1)
%         col_Hi=H(:,j);
%         Num_nonzero_H(j)=size(find(col_Hi==T1_hop|col_Hi==T2_hop),1);
%     end
    %Num_nonzero_H
%     nonzeros_H=zeros(max(Num_nonzero_H),size(H,1));
%     for j=1:size(H,1)
%         col_Hi=H(:,j);
%         nonzero_H_i=find(col_Hi==T1_hop|col_Hi==T2_hop);
%         nonzeros_H(:,j)=[nonzero_H_i;zeros(max(Num_nonzero_H)-size(nonzero_H_i,1),1)];

%     end
%     [1:size(H,1);nonzeros_H]
    
     Num_nonzero_H_spin1=zeros(1,size(H_spin1,1));
     for j=1:size(H_spin1,1)
         col_Hi=H_spin1(:,j);
         Num_nonzero_H_spin1(j)=size(find(col_Hi==T1_hop),1);
     end
    %Num_nonzero_H_spin1
     nonzeros_H_spin1=zeros(max(Num_nonzero_H_spin1),size(H_spin1,1));
     for j=1:size(H_spin1,1)
         col_Hi=H_spin1(:,j);
         nonzero_H1_i=find(col_Hi==T1_hop);
         nonzeros_H_spin1(:,j)=[nonzero_H1_i;zeros(max(Num_nonzero_H_spin1)-size(nonzero_H1_i,1),1)];
     end
     
     Hspin1_bond=[1:size(H_spin1,1);nonzeros_H_spin1];
     hbonds_translation_spin1=HBond_translation(bases_fourier_sector,Hspin1_bond);
     for khbonds=1:size(bases_fourier_sector.num,2)
         hbonds_translation_spin1{khbonds}
     end
    
     Num_nonzero_H_spin2=zeros(1,size(H_spin2,1));
     for j=1:size(H_spin2,1)
         col_Hi=H_spin2(:,j);
         Num_nonzero_H_spin2(j)=size(find(col_Hi==T2_hop),1);
     end
    %Num_nonzero_H_spin2
     nonzeros_H_spin2=zeros(max(Num_nonzero_H_spin2),size(H_spin2,1));
     for j=1:size(H_spin2,1)
         col_Hi=H_spin2(:,j);
         nonzero_H2_i=find(col_Hi==T2_hop);
         nonzeros_H_spin2(:,j)=[nonzero_H2_i;zeros(max(Num_nonzero_H_spin2)-size(nonzero_H2_i,1),1)];
     end
     
     Hspin2_bond=[1:size(H_spin2,1);nonzeros_H_spin2];
     hbonds_translation_spin2=HBond_translation(bases_fourier_sector,Hspin2_bond);
     for khbonds=1:size(bases_fourier_sector.num,2)
         hbonds_translation_spin2{khbonds}
     end
    
    momentum_max = Nsites-1;
    for momentum=0:momentum_max
        momentum
        if size(bases_fourier_sector.num(bases_fourier_sector.num<Nsites),2)>0
            states_momentum=States_momentum(bases_fourier_sector,momentum);
        else
            states_momentum=1:size(bases_fourier_sector.num,2);
        end
        states_momentum
        H_fourier_spin1=Hamiltonian_fourier(momentum,states_momentum,hbonds_translation_spin1,bases_fourier_sector,1);
        H_fourier_spin2=Hamiltonian_fourier(momentum,states_momentum,hbonds_translation_spin2,bases_fourier_sector,2);
        H_fourier=H_fourier_spin1+H_fourier_spin2;
        %[eigenstates_f,energies]=eig(H_fourier);
        %energy_list=diag(energies)
        %eigenstates_f
        energy_list=eig(H_fourier);
        e_diff=diff(energy_list);
        e_diff=e_diff./mean(e_diff);
%         disp("Number of zeros: ")
%         disp(sum(e_diff==0))
    
%         len_e_diff=size(e_diff,1);
%         r_ediff=(e_diff(2:len_e_diff))./(e_diff(1:len_e_diff-1));
%        binwidth=min(nonzeros(r_ediff));
%         binwidth=min(nonzeros(e_diff));
        figure(fignum)
        fignum=fignum+1;
        %hi=histogram(energy_list);
%         hi=histogram(r_ediff,'Normalization','pdf');
        hi=histogram(e_diff,'Normalization','pdf');
%         hi.BinWidth=binwidth;
        xmax_plot=5;
        hi.BinLimits=[0,xmax_plot];
        x_plot=0:0.05:xmax_plot;
        hold on 
        plot(x_plot,poisson(x_plot),'-.r')
        hold on
        plot(x_plot,GOE(x_plot),'--b')
        title(['level stat, sector: ',num2str(sector), '(', num2str([Nsites, Nparticles1, Nparticles2]),'),', ' momentum= ', num2str(momentum), '(2\pi/',num2str(Nsites),'a)'])
        saveas(hi,[pwd '/plots/1477/without inv sym/level stat/' num2str(sector) '_p=' num2str(momentum) '.pdf']);
    end
    
end
end



function f=poisson(r)
%Poisson distribution of r_n=s_n/s_n-1
%f=1./((1+r).^2);
f=exp(-r);
end

function f=GOE(r)
% Wigner-Dyson distribution GOE 
% Z=8/27;
% f=(1/Z).*(r+r.^2)./((1+r+r.^2).^2.5);

f=(pi./2).*r.*exp(-(pi./4).*r.^2);
end

function hbonds_translation=HBond_translation(bases_f,H_bond)
%bases_f is of class BasesFourier
%only get the upper triangular elements
dim=size(bases_f.num,2);
hbonds_translation=cell(1,dim);
for ind_states=1:dim
     bondsf_i=Bonds_transym;
     bondsf_i.stateind=ind_states;
     bondsf_i.bonds=[];
     bondsf_i.offsets=[];
     statei_f=bases_f.states{ind_states};
     statei=statei_f(1);
     bondi=H_bond(:,statei);
     bondi=transpose(bondi(bondi>statei));
     if size(bondi,1)~=0 && size(bondi,2)~=0
         for statej=bondi
             [ind_statejf, offsetj]=find_state_in_basesfourier(statej,bases_f);
             if ~ismember(ind_statejf,bondsf_i.bonds)
                 bondsf_i.bonds=[bondsf_i.bonds,ind_statejf];
                 bondsf_i.offsets=[bondsf_i.offsets,offsetj];
             end
         end
     end
     hbonds_translation{ind_states}=bondsf_i;
    
end
end


function [ind_statef, offset]=find_state_in_basesfourier(state,bases_f)
%bases_f is of class BasesFourier

dimh=size(bases_f.num,2);
firsts=zeros(1,dimh);
for i=1:dimh
    firsts(i)=bases_f.states{i}(1);
end
largergroup=find(firsts>state);
if size(largergroup,1)==0 || size(largergroup,2)==0 
    ind_statef=dimh;
else
    ind_statef=largergroup(1)-1;
end
offset=state-firsts(ind_statef);
end



function H_fourier=Hamiltonian_fourier(momentum,states_momentum,hbonds_translation,bases_f,spin)
%hbonds_translation is of class Bonds_transym
%matrix elements of two vector states l1/l2==m is
%sqrt(m).T.(omega^(momentum))^offset
%spin==1 or spin==2
global omega T1_hop T2_hop
unitel=omega^momentum;
dimh=size(states_momentum,2);
H_fourier=zeros(dimh,dimh);
for i=1:dimh
    statei_f=states_momentum(i);
    leni=bases_f.num(statei_f);
    
    [statej_f_choices,ia]=intersect(hbonds_translation{statei_f}.bonds,states_momentum);
    offsetsi_f_choices= hbonds_translation{statei_f}.offsets(ia);
    num_j_choices=size(statej_f_choices,2);
    if ~ismember(0,size(statej_f_choices))
    for jj=1:num_j_choices
        statej_f=statej_f_choices(jj);
        j=find(states_momentum==statej_f);
        lenj=bases_f.num(statej_f);
        if leni==lenj
            H_fourier(i,j)=unitel^offsetsi_f_choices(jj);
        else
            l1=max(leni,lenj);
            l2=min(leni,lenj);
            if mod(l1,l2)~=0
                disp('error')
                pause
            else
                m=l1/l2;
                H_fourier(i,j)=sqrt(m)*(unitel^offsetsi_f_choices(jj));
            end
        end
    end
    end
end
H_fourier=H_fourier'+H_fourier;
if spin==1
    H_fourier=H_fourier*T1_hop;
else
    H_fourier=H_fourier*T2_hop;
end
end



function states=States_momentum(bases_f,momentum)
states=[];
for i=1:size(bases_f.num,2)
    if ismember(momentum,bases_f.momenta{i})
        states=[states,i];
    end
end
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








function [stateMatrix,bases_fourier]=States_sector(indspin1_start, sector)
% sum(sector)=Nparticles2
% length(sector)=Nparticles1
% every two rows represent one state (bot row--spin 1;top row--spin 2)
global Nsites Nparticles1 


%for a given sector, how choices for the indices of spin 1
len=Nparticles1; %len is the number of boxes
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

indspin1 = DeleteIndspin1Duplicates(sector,indspin1_seed); % each row of indspin1 tells the sizes of #len boxes
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
bases_fourier=BasesFourier;
bases_fourier.num=[];
bases_fourier.states={};
num_stateMatrix=0;
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
    
   
    if maxSteps_circshift(i)==Nsites
        for num_spin2_i=1:size(states_spin2_i,1)
            state_combined_1_num=[states_spin2_i(num_spin2_i,:);states_spin1(i,:)];
            bases_fourier.num=[bases_fourier.num,Nsites];
            basis_f=zeros(1,Nsites);
            for shiftstep=0:Nsites-1
                stateMatrix=[stateMatrix;circshift(state_combined_1_num,shiftstep,2)];
                num_stateMatrix=num_stateMatrix+1;
                basis_f(shiftstep+1)=num_stateMatrix;
            end
            bases_fourier.states=[bases_fourier.states,basis_f];
        end
    else
        pool=1:size(states_spin2_i,1);
        selected={};%{[period,state_start]}
        while size(pool,2)>0
            num_p1=pool(1);
            state_pair1=[states_spin2_i(num_p1,:);states_spin1(i,:)];
            
            selected_p1=[0,num_p1];
            
            shiftstep_unit=maxSteps_circshift(i);
            shiftstep_max=Nsites/shiftstep_unit;
            for shiftstep_m=1:shiftstep_max
                shiftstep=shiftstep_unit*shiftstep_m;
                state_pair1_shift=circshift(state_pair1,shiftstep,2);
                for num_p2=pool
                    state_pair2=[states_spin2_i(num_p2,:);states_spin1(i,:)];
                    if sum(sum(state_pair1_shift==state_pair2))==2*Nsites
                        pool(pool==num_p2)=[];
                        break
                    end
                end
                if num_p2==num_p1
                    selected_p1(1)=shiftstep;
                    break
                end
            end
            selected=[selected,selected_p1];
        end
        
        for num_spin2_i=1:size(selected,2)
            state_combined_1_num=[states_spin2_i(selected{num_spin2_i}(2),:);states_spin1(i,:)];
            period=selected{num_spin2_i}(1);
            bases_fourier.num=[bases_fourier.num,period];
            basis_f=zeros(1,period);
            for shiftstep=0:period-1
                stateMatrix=[stateMatrix;circshift(state_combined_1_num,shiftstep,2)];
                num_stateMatrix=num_stateMatrix+1;
                basis_f(shiftstep+1)=num_stateMatrix;
            end
            bases_fourier.states=[bases_fourier.states,basis_f];
        end 
    end
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
