function relax_dynamics_energy1
%To compare the relaxation times of 1D Fermi gas and Haldane (semion) particle gas
%with linear dispersion relation \epsilon=v*k
%coupled a bath of phonons
%energy levels are equally spaced and spacing between nearest levels is 1
%Lindblad operator:
%L1--hopping from high-energy level to low-energy level,
%L2--hopping from low-energy level to high-energy level.
%Hermitian conjugate of L1 is L2


clear all

Nlevel=6;
Nparticle_fspinless=3;
Nparticle=2*Nparticle_fspinless;
fermilevel=Nparticle_fspinless;

global T gamma imagunit dt
imagunit=sqrt(-1);
T=0.01; %temperature
gamma=10;%gamma=2\pi*Density of states*(coupling constant)^2
dt=0.01;


nmax=200;
tvector=0:dt:(nmax-1)*dt;

stateMatrixFermi=FermiStates(Nlevel,Nparticle);
[row_f,col_f]=size(stateMatrixFermi);
stateMatrixSemion=SemionStates(Nlevel, Nparticle)
[row_s,col_s]=size(stateMatrixSemion)
stateMatrixFermi_spinless=FermiStates_spinless(Nlevel,Nparticle_fspinless)
[row_fspl,col_fspl]=size(stateMatrixFermi_spinless);


%levels to be considered
levels_consider=[fermilevel-1,fermilevel,fermilevel+1];
sizelevel_consider=size(levels_consider);
clevel_consider=sizelevel_consider(2);

%particle move by nemax level units
nemax=2;



%H0 of spinful fermion
% H0_Fermion_spinful=H0_spin(stateMatrixFermi,row_f,col_f);

%H0 of spinless fermion
H0_Fermion_spinless=H0_spinless(stateMatrixFermi_spinless,row_fspl,col_fspl);

%H0 of semion
H0_Semion=H0_spin(stateMatrixSemion,row_s,col_s);






%matrix form of occpuation number operator of spinful fermion
% occupationMatrixFermi=zeros(row_f,row_f,clevel_consider);
% for i=1:clevel_consider
%     levelis=levels_consider(i);
%     nmatrixis=OccupationMatrix_semion(stateMatrixFermi,row_f,levelis);
%     occupationMatrixFermi(:,:,i)=nmatrixis;
% end


%matrix form of occpuation number operator of spinless Fermion 
occupationMatrixFermi_spinless=zeros(row_fspl,row_fspl,clevel_consider);
for i=1:clevel_consider
    leveli=levels_consider(i);
    nmatrixi=OccupationMatrix_spinlessfermion(stateMatrixFermi_spinless,row_fspl,leveli);
    occupationMatrixFermi_spinless(:,:,i)=nmatrixi;
end
%occupationMatrixFermi_spinless


%matrix form of occpuation number operator of Semion
occupationMatrixSemion=zeros(row_s,row_s,clevel_consider);
for i=1:clevel_consider
    levelis=levels_consider(i);
    nmatrixis=OccupationMatrix_semion(stateMatrixSemion,row_s,levelis);
    occupationMatrixSemion(:,:,i)=nmatrixis;
end
%occupationMatrixSemion





%Gammaton
Gamma_ton=zeros(nemax,2);
for i=1:nemax
    Gamma_ton(i,:)=Gamma(i);
end


%Lindblad operator of spinful Fermion
% L2ton_FermionSpinful=zeros(row_f,row_f,nemax);
% for i=1:nemax
%     L2fermionSpinful_i=Lindblad_Fermi(stateMatrixFermi,row_f,i);
%     L2ton_FermionSpinful(:,:,i)=L2fermionSpinful_i;
% end
% L2ton_FermionSpinful(:,:,2)

%Lindblad operator of spinless Fermion
L2ton_Fermion=zeros(row_fspl,row_fspl,nemax);
for i=1:nemax
    L2fermion_i=Lindblad_Fermi(stateMatrixFermi_spinless,row_fspl,i);
    L2ton_Fermion(:,:,i)=L2fermion_i;
end
% L2ton_Fermion(:,:,1)

%Lindblad operator of Semion
L2ton_Semion=zeros(row_s,row_s,nemax);
for i=1:nemax
    L2semion_i=Lindblad_Semion(stateMatrixSemion, row_s, i);
    L2ton_Semion(:,:,i)=L2semion_i;
end
% L2ton_Semion(:,:,1)



%Initial state/density matrix of spinful Fermion
% Rho_f_spinful=zeros(row_f,row_f);
% Rho_f_spinful(row_f,row_f)=1;


%Initial state/density matrix of spinless Fermion
Rho_f=zeros(row_fspl,row_fspl);
%nrho_fspl=0;
%for i=1:row_fspl
%    vfi=stateMatrixFermi_spinless(i,:);
%    if IdenticalVectors(stateInitial_s,vfi,col_fspl)==1


%        nrho_f=i;
%    end
%end
%Rho_f(nrho_f,nrho_f)=1;
%Rho_f(row_fspl-6,row_fspl-6)=1;
%Rho_f(1,1)=1;
%Rho_f(3,6)=imagunit;
%Rho_f(6,3)=-imagunit;
for kk=1:row_fspl
    Rho_f(kk,kk)=random('Uniform',0,1);
end
trace_Rho_f=trace(Rho_f);
for kk=1:row_fspl
    Rho_f(kk,kk)=Rho_f(kk,kk)./trace_Rho_f;
end


for kr=1:row_fspl
    for kc=kr+1:row_fspl
        re = random('Uniform',0,1);
        imag = random('Uniform',0,1);
        Rho_f(kr,kc) = re + imagunit.*imag;
        Rho_f(kc,kr) = re - imagunit.*imag;
    end
end


%Initial state/density matrix of Semion
Rho_s=zeros(row_s,row_s);

for kk=1:row_s
    Rho_s(kk,kk)=random('Uniform',0,1);
end
trace_Rho_s=trace(Rho_s);
for kk=1:row_s
    Rho_s(kk,kk)=Rho_s(kk,kk)./trace_Rho_s;
end

for kr=1:row_s
    for kc=kr+1:row_s
        re = random('Uniform',0,1);
        imag = random('Uniform',0,1);
        Rho_s(kr,kc) = re + imagunit.*imag;
        Rho_s(kc,kr) = re - imagunit.*imag;
    end
end


%Rho_s(row_s-17,row_s-17)=1;
%Rho_s(1,1)=1;
%Rho_s(3,6)=imagunit;
%Rho_s(6,3)=-imagunit;

%Initial occpuation number of spinful Fermion
% occupationNumberFermi_spinful=zeros(clevel_consider,nmax);
% for k=1:clevel_consider
%     initialocci=trace(Rho_f_spinful*occupationMatrixFermi(:,:,k));
%     occupationNumberFermi_spinful(k,1)=initialocci;    
% end



%Initial occpuation number of spinless Fermion
occupationNumberFermi=zeros(clevel_consider,nmax);
for k=1:clevel_consider
    initialocci=trace(Rho_f*occupationMatrixFermi_spinless(:,:,k));
    occupationNumberFermi(k,1)=initialocci;    
end


%Initial occpuation number of Semion
occupationNumberSemion=zeros(clevel_consider,nmax);
for k=1:clevel_consider
    initialocci_s=trace(Rho_s*occupationMatrixSemion(:,:,k));
    occupationNumberSemion(k,1)=initialocci_s;    
end


%Energy of  initial state of spinful fermion
% energyVectorFermi_spinful=zeros(1,nmax);
% energyVectorFermi_spinful(1)=trace(Rho_f_spinful*H0_Fermion_spinful);


%Energy of  initial state of spinless fermion
energyVectorFermi=zeros(1,nmax);
energyVectorFermi(1)=trace(Rho_f*H0_Fermion_spinless);


%Energy of  initial state of semion
energyVectorSemion=zeros(1,nmax);
energyVectorSemion(1)=trace(Rho_s*H0_Semion);


densityMatrix_f_diag=zeros(row_fspl,nmax);
densityMatrix_f_off_diag=zeros(row_fspl.*(row_fspl-1)./2,nmax);

for kk=1:row_fspl
    densityMatrix_f_diag(kk,1)=Rho_f(kk,kk);
end

k_count = 0;
for kr=1:row_fspl
    for kc=kr+1:row_fspl
        k_count = k_count + 1;
        densityMatrix_f_off_diag(k_count,1)=Rho_f(kr,kc);
    end
end


densityMatrix_s_diag=zeros(row_s,nmax);
densityMatrix_s_off_diag=zeros(row_s.*(row_s-1)./2,nmax);

for kk=1:row_s
    densityMatrix_s_diag(kk,1)=Rho_s(kk,kk);
end

k_count = 0;
for kr=1:row_s
    for kc=kr+1:row_s
        k_count = k_count + 1;
        densityMatrix_s_off_diag(k_count,1)=Rho_s(kr,kc);
    end
end

 %Rho_s36=zeros(1,nmax);
 %Rho_s33=zeros(1,nmax);
 %Rho_f36=zeros(1,nmax);
 %Rho_f33=zeros(1,nmax);
% 
 %Rho_s36(1)=Rho_s(3,6);
 %Rho_s33(1)=Rho_s(3,3);
 %Rho_f36(1)=Rho_f(3,6);
 %Rho_f33(1)=Rho_f(3,3);

 
%Calculate next denstity matrix using Lindblad equation

for i=2:nmax
    
%    Rho_f_spinful=NextDensityMatrix(H0_Fermion_spinful,L2ton_FermionSpinful,Gamma_ton,Rho_f_spinful,nemax);
   Rho_f=NextDensityMatrix(H0_Fermion_spinless,L2ton_Fermion,Gamma_ton,Rho_f,nemax);
   Rho_s=NextDensityMatrix(H0_Semion,L2ton_Semion,Gamma_ton,Rho_s,nemax);
   
   
   for kk=1:row_fspl
    densityMatrix_f_diag(kk,i)=Rho_f(kk,kk);
   end
   
   k_count = 0;
   for kr=1:row_fspl
       for kc=kr+1:row_fspl
           k_count = k_count + 1;
           densityMatrix_f_off_diag(k_count,i)=Rho_f(kr,kc);
       end
   end

   for kk=1:row_s
       densityMatrix_s_diag(kk,i)=Rho_s(kk,kk);
   end
   
   k_count = 0;
   for kr=1:row_s
       for kc=kr+1:row_s
           k_count = k_count + 1;
           densityMatrix_s_off_diag(k_count,i)=Rho_s(kr,kc);
       end
   end
   
   
    %Rho_s36(i)=Rho_s(3,6);
    %Rho_s33(i)=Rho_s(3,3);
    %Rho_f36(i)=Rho_f(3,6);
    %Rho_f33(i)=Rho_f(3,3);
   
   for k=1:clevel_consider
%        occupationNumberFermi_spinful(k,i)=trace(Rho_f_spinful*occupationMatrixFermi(:,:,k));
       occupationNumberFermi(k,i)=trace(Rho_f*occupationMatrixFermi_spinless(:,:,k));
       occupationNumberSemion(k,i)=trace(Rho_s*occupationMatrixSemion(:,:,k));
   end
   
%    energyVectorFermi_spinful(i)=trace(Rho_f_spinful*H0_Fermion_spinful);
   energyVectorFermi(i)=trace(Rho_f*H0_Fermion_spinless);
   energyVectorSemion(i)=trace(Rho_s*H0_Semion);
end



%plots of density matrix elements vs time
figure_num=1;
figure(figure_num)
figure_num = figure_num + 1;
for kk=1:row_fspl
    %figure(figure_num);
    %figure_num = figure_num+1;
    
    plot(tvector, densityMatrix_s_diag(kk,:),'b');
    hold on
end
str = ["Semion diagonal elements"];
title(str);
%legend('Fermion')

figure(figure_num);
%figure_num = figure_num + 1;
k_count = 0;
for kr=1:row_fspl
    for kc=kr+1:row_fspl
        k_count = k_count + 1;
        plot(tvector, abs(densityMatrix_s_off_diag(k_count,:)),'b');
        hold on
    end
end
str = ["Semion off diagonal elements"];
title(str);
%legend('Fermion')

%plots of occupation number vs time

% for k=1:3
%     figure(k)
%     plot(tvector,2*occupationNumberFermi(k,:),'r-.')
%     hold on
%     plot(tvector,occupationNumberSemion(k,:),'b')
%     str=['Fermi level (+) ',num2str(levels_consider(k)-fermilevel)];
%     title(str)
%     hold off
%     legend('Spinful fermion','Semion')
% end




% semionEXP=[2.,1.99998,0.000025];
% fermionEXP=[2.,1.98661,0.0133855];
% for k=1:2
% %     figure(k)
% %     plot(tvector,2*occupationNumberFermi(k,:),'r-.')
% %     hold on
% %     plot(tvector,occupationNumberSemion(k,:),'b')
% %     hold on
% %     plot(tvector,occupationNumberFermi_spinful(k,:),'g--')
% %     str=['Fermi level (+) ',num2str(levels_consider(k)-fermilevel)];
% %     title(str)
% %     hold off
% %     legend('Spinful fermion','Semion')
% %     legend('Spinless fermion(double)','Semion','Spinful fermion')
%     
%     
%     
%     
%     figure(k)
%     plot(tvector,(fermionEXP(k)/(2*occupationNumberFermi(k,nmax)))*2*occupationNumberFermi(k,:),'r-.')
%     hold on
%     plot(tvector,(semionEXP(k)/occupationNumberSemion(k,nmax))*occupationNumberSemion(k,:),'b')
   
    
%     semionEXPk=semionEXP(k)*ones(1,nmax);
%     fermionEXPk=fermionEXP(k)*ones(1,nmax);
%     hold on
%     plot(tvector,fermionEXPk,'r-.')
%     hold on
%     plot(tvector,semionEXPk,'b-.')

%     str=['Fermi level (+) ',num2str(levels_consider(k)-fermilevel)];
%     title(str)
%     hold off
%     legend('Spinful fermion','Semion')
% end

% k=k+1;
% nmid_f=0.045/dt;
% max_f=max(2*occupationNumberFermi(k,1:nmax));
% exp_f=fermionEXP(3);
% equilib_f3=2*occupationNumberFermi(k,nmax);
% occupationNumberFermi1st=2*occupationNumberFermi(k,1:nmid_f);
% occupationNumberFermi2nd=2*occupationNumberFermi(k,nmid_f+1:nmax);
% occupationNumberFermi2nd=(occupationNumberFermi2nd-max_f)*((max_f-exp_f)/(max_f-equilib_f3))+max_f;
% occupationNumberFermicorr=[occupationNumberFermi1st,occupationNumberFermi2nd];
% nmid_s=0.415/dt;
% max_s=max(occupationNumberSemion(k,1:nmax));
% exp_s=semionEXP(3);
% equilib_s3=occupationNumberSemion(k,nmax);
% occupationNumberSemion1st=occupationNumberSemion(k,1:nmid_s);
% occupationNumberSemion2nd=occupationNumberSemion(k,nmid_s+1:nmax);
% occupationNumberSemion2nd=(occupationNumberSemion2nd-max_s)*((max_s-exp_s)/(max_s-equilib_s3))+max_s;
% occupationNumberSemioncorr=[occupationNumberSemion1st,occupationNumberSemion2nd];
% figure(k)
% plot(tvector,occupationNumberFermicorr,'r-.')
% hold on
% plot(tvector,occupationNumberSemioncorr,'b')
% str=['Fermi level (+) ',num2str(levels_consider(k)-fermilevel)];
% title(str)
% hold off
% legend('Spinful fermion','Semion')



%plot of energy vs time

% figure(k+1)
% plot(tvector,2*energyVectorFermi,'r-.')
% hold on
% plot(tvector,energyVectorSemion,'b')
% % hold on
% % plot(tvector,energyVectorFermi_spinful,'g--')
% title("energy vs time")
% hold off
% legend('Spinful fermion','Semion')
% % legend('Spinless fermion(double)','Semion','Spinful fermion')

% energyG=6;
% energyEXP_f=energyG+0.0133841;
% energyEXP_s=energyG+0.000025;
% energyFermi=(2*energyVectorFermi-2*energyVectorFermi(1))*((2*energyVectorFermi(1)-energyEXP_f)/(2*energyVectorFermi(1)-2*energyVectorFermi(nmax)))+2*energyVectorFermi(1);
% energySemion=(energyVectorSemion-energyVectorSemion(1))*((energyVectorSemion(1)-energyEXP_s)/(energyVectorSemion(1)-energyVectorSemion(nmax)))+energyVectorSemion(1);
% figure(k+1)
% plot(tvector,energyFermi,'r-.')
% hold on
% plot(tvector,energySemion,'b')
% title("energy vs time")
% hold off
% legend('Spinful fermion','Semion')


%plot of log(energy) vs time

% endp=nmax;
% figure(k+2)
% plot(tvector(1:endp),log(2*energyVectorFermi(1:endp)),'r-.')
% hold on
% plot(tvector(1:endp),log(energyVectorSemion(1:endp)),'b')
% title("log(energy) vs time")
% legend('Fermion','Semion')
% hold off
% 
% 
% endp=nmax/2;
% figure(k+3)
% plot(tvector(1:endp),log(2*energyVectorFermi(1:endp)),'r-.')
% hold on
% plot(tvector(1:endp),log(energyVectorSemion(1:endp)),'b')
% title("log(energy) vs time")
% legend('Fermion','Semion')
% hold off



% figure(k+6)
% plot(tvector(1:endp-1),abs(log(2*energyVectorFermi(2:endp))-log(2*energyVectorFermi(1:endp-1)))./dt,'r-.')
% hold on
% plot(tvector(1:endp-1),abs(log(energyVectorSemion(2:endp))-log(energyVectorSemion(1:endp-1)))./dt,'b')
% title("dlog(energy)/dt vs time")
% legend('Fermion','Semion')
% hold off

% figure(k+7)
% plot(energyFermi(2:endp),abs(log(2*energyVectorFermi(2:endp))-log(2*energyVectorFermi(1:endp-1)))./dt,'r-.')
% hold on
% plot(energySemion(2:endp),abs(log(energyVectorSemion(2:endp))-log(energyVectorSemion(1:endp-1)))./dt,'b')
% title("dlog(energy)/dt vs energy")
% legend('Spinful fermion','Semion')
% hold off

% figure(k+7)
% plot(2*energyVectorFermi(1:endp-1)-6,abs(log(2*energyVectorFermi(2:endp))-log(2*energyVectorFermi(1:endp-1)))./dt,'r-.')
% 
% hold on
% plot(energyVectorSemion(1:endp-1)-6,abs(log(energyVectorSemion(2:endp))-log(energyVectorSemion(1:endp-1)))./dt,'b')
% title("dlog(energy)/dt vs energy")
% legend('Spinful fermion','Semion')
% hold off

% hold on
% plot((24/(max(energyVectorFermi_spinful)-min(energyVectorFermi_spinful)))*(energyVectorFermi_spinful(1:endp-1)-energyVectorFermi_spinful(endp-1))+6,abs(log(energyVectorFermi_spinful(2:endp))-log(energyVectorFermi_spinful(1:endp-1)))./dt,'g--')
% legend('Spinless fermion(double)','Semion','Spinful fermion')



%  figure(k+8)
%  plot(tvector,real(Rho_s36),'r-.')
%  hold on
%  plot(tvector,imag(Rho_s36),'g--')
%  hold on
%  plot(tvector,Rho_s33,'b')
%  title("density matrix elements (Haldane)")
%  legend('real_Rho_s36','imag_Rho_s36','Rho_s33')
%  hold off
% % 
% % 
%  figure(k+9)
%  plot(tvector,real(Rho_f36),'r-.')
%  hold on
%  plot(tvector,imag(Rho_f36),'g--')
%  hold on
%  plot(tvector,Rho_f33,'b')
%  title("density matrix elements (Fermion)")
%  legend('real_Rho_f36','imag_Rho_f36','Rho_f33')
%  hold off




end



function stateMatrix=FermiStates_spinless(nlevel,nparticle_spinless)
v=1:1:nlevel;
state_spinless=nchoosek(v,nparticle_spinless);
[nrow,ncolumn]=size(state_spinless);
stateMatrix=zeros(nrow,nlevel);
for i=1:nrow
    for j=1:ncolumn
        m=state_spinless(i,j);
        stateMatrix(i,m)=1;
    end
end

end




function stateMatrix=FermiStates(nlevel, nparticle)
% states are of the form (0 up, 0 down, 1 up, 1 down,..., nlevel up, nlevel down)

v=1:1:nlevel;
n=nparticle./2;
state_spinless=nchoosek(v,n);
[nrow,ncolumn]=size(state_spinless);
nrow_state=nrow*nrow;
stateMatrix=zeros(nrow_state,2.*nlevel);
for i=1:nrow
    nrs=(i-1)*nrow;
    for j=1:nrow
        nrrs=nrs+j;
        
        for k1=1:ncolumn
            m1=2*state_spinless(i,k1)-1;
            stateMatrix(nrrs,m1)=1;
        end
        for k2=1:ncolumn
            m2=2*state_spinless(j,k2);
            stateMatrix(nrrs,m2)=1; 
        end
    end
end

end



function stateMatrix=SemionStates(nlevel, nparticle)
% To get rid of Fermi states which don't obey the occupation rules of
% Semion Haldane statistics (g=1/2)

v=1:1:nlevel;
n=nparticle./2;
state_spinless=nchoosek(v,n);
[nrow,ncolumn]=size(state_spinless);

stateMatrix=[];
for i=1:nrow
    for j=1:nrow
        %v1 is the lower row, or flavor 1; v2 is the upper row, or flavor 2
        v1=state_spinless(i,:);
        v2=state_spinless(j,:);
        flag=SemionFilter(v1,v2,ncolumn,nlevel);
        
        if flag==1
            vs=zeros(1,2*nlevel);
            
            for k1=1:ncolumn
                m1=2*v1(k1)-1;
                vs(m1)=1;
            end
            
            for k2=1:ncolumn
                m2=2*v2(k2);
                vs(m2)=1;
            end
            
            stateMatrix=[stateMatrix;vs];
            
        else
            continue
        end
    end
end

end


function flag=SemionFilter(v1,v2,vsize,maxlevel)
%to filter out states which don't obey Semion occupation constraints
v1=[v1,maxlevel+1];
flag=1;
for i=1:vsize
    if v2(i)<v1(i)||v2(i)>=v1(i+1)
        flag=0;
        break
    end
end
end




function Hmatrix=H0_spin(stateMatrix,row,col)
%H0 of many-body states for semion or spinful fermion
Hmatrix=zeros(row,row);
energyLevels=zeros(1,col);
for i=1:(col/2)
    energyLevels((i-1)*2+1)=i-1;
    energyLevels((i-1)*2+2)=i-1;
end

for i=1:row
    vi=stateMatrix(i,:);
    ei=sum(vi.*energyLevels);
    Hmatrix(i,i)=ei;
end
end


function Hmatrix=H0_spinless(stateMatrix,row,col)
%H0 of many-body states for spinless fermion
Hmatrix=zeros(row,row);
energyLevels=0:1:col-1;
for i=1:row
    vi=stateMatrix(i,:);
    ei=sum(vi.*energyLevels);
    Hmatrix(i,i)=ei;
end
end


function nmatrix=OccupationMatrix_semion(stateMatrix,row,n)
%Matrix of occupation number of level n for semion
nmatrix=zeros(row,row);
for i=1:row
    vi=stateMatrix(i,:);
    nmatrix(i,i)=vi((n-1)*2+1) + vi((n-1)*2+2);
end
end


function nmatrix=OccupationMatrix_spinlessfermion(stateMatrix,row,n)
%Matrix of occupation number of level n for fermion
nmatrix=zeros(row,row);
for i=1:row
    vi=stateMatrix(i,:);
    nmatrix(i,i)=vi(n);    
end
end


function n_b=BoseFun(n)
%n:w=n*energy units
global T
n_b=1./(exp(n./T)-1);
end


function Gamma_n=Gamma(n)
global gamma
n_b=BoseFun(n);
Gamma1=gamma.*n_b;
Gamma2=gamma.*(n_b+1);
Gamma_n=[Gamma1,Gamma2];
end


function Ln_Fermi=Lindblad_Fermi(stateMatrixFermi,row,n)
%Lindblad operator for Fermions moving by n energy units
%hop from high-energy level to low-energy level
Ln_Fermi=zeros(row,row);
for i=1:row-1
    for j=i+1:row
        vi=stateMatrixFermi(i,:);
        vj=stateMatrixFermi(j,:);
        [flag,p_particle,p_hole]=hopFromjToi_fermi(vi,vj,n);
        if flag==1
            vbetween=vj(p_hole+1:p_particle-1);
            Lij=(-1)^(sum(vbetween));
            Ln_Fermi(i,j)=Lij;
        elseif flag==-1
            vbetween=vi(p_hole+1:p_particle-1);
            Lji=(-1)^(sum(vbetween));
            Ln_Fermi(j,i)=Lji;
        end
        %pause
    end
end
end


function [flag,p_particle,p_hole]=hopFromjToi_fermi(vi,vj,n)
%hop from high-energy level to low-energy level
%hopping between levels with energy difference being n 
%flag=0 no hopping can happen between state i and j
%flag=1 can hop from state j to i
%flag=-1 can hop from state i to j

flag=0;
p_particle=0;
p_hole=0;
deltav=vj-vi;
if sum(abs(deltav))==2
    p_particle=find(deltav==1);
    p_hole=find(deltav==-1);
    deltap=p_particle-p_hole;
    if abs(deltap)==n
        if deltap>0
            flag=1;
        else
            flag=-1;
        end
    end 
end
end



function Ln_Semion=Lindblad_Semion(stateMatrixSemion, row, n)
%Lindblad operator for Semions moving by n energy units
%hop from high-energy level to low-energy level
Ln_Semion=zeros(row,row);
for i=1:row-1
    for j=i+1:row
        vi=stateMatrixSemion(i,:);
        vj=stateMatrixSemion(j,:);
        flag=hopFromjToi_semion(vi,vj,n);
        if flag==1
            Ln_Semion(i,j)=1;
        elseif flag==-1
            Ln_Semion(j,i)=1;
        end
        %pause
    end
end


end




function flag=hopFromjToi_semion(vi,vj,n)
%hop from high-energy level to low-energy level
%hopping between levels with energy difference being n 
%flag=0 no hopping can happen between state i and j
%flag=1 can hop from state j to i
%flag=-1 can hop from state i to j

flag=0;
deltav=vj-vi;
if sum(abs(deltav))==2
    p_particle=find(deltav==1);
    p_hole=find(deltav==-1);
    deltap=p_particle-p_hole;
    if (abs(deltap)/2)==n
        if deltap>0
            flag=1;
        else
            flag=-1;
        end
    end 
end
end




function densityMatrix_next=NextDensityMatrix(H0matrix,L2ton,Gamma_ton,densityMatrix,n)
%To find density matrix of next step in Runge-Kutta algorithm
global dt
K1=LindbladEqn(H0matrix,L2ton,Gamma_ton,densityMatrix,n);
K2=LindbladEqn(H0matrix,L2ton,Gamma_ton,densityMatrix+0.5*dt*K1,n);
K3=LindbladEqn(H0matrix,L2ton,Gamma_ton,densityMatrix+0.5*dt*K2,n);
K4=LindbladEqn(H0matrix,L2ton,Gamma_ton,densityMatrix+dt*K3,n);
densityMatrix_next=densityMatrix + (dt/6)*(K1+2*K2+2*K3+K4);
end



function fmatrix=LindbladEqn(H0matrix,L2ton,Gamma_ton,densityMatrix,n)
%make the right-hand side of Lindblad equation
global imagunit
dimh = size(densityMatrix,1);
fmatrix=-imagunit*Commutator(H0matrix,densityMatrix);
for i=1:n
    Gamma_i=Gamma_ton(i,:);
    L2_i=L2ton(:,:,i);
    for j=1:dimh-1
        for k=j+1:dimh
            if L2_i(j,k)~=0
                L2_is = zeros(dimh,dimh);
                L2_is(j,k) = L2_i(j,k);
                L1_is=transpose(L2_is);
                dissipation1 = Gamma_i(1)*(L1_is*densityMatrix*transpose(L1_is)-0.5*antiCommutator(transpose(L1_is)*L1_is,densityMatrix));
                dissipation2 = Gamma_i(2)*(L2_is*densityMatrix*transpose(L2_is)-0.5*antiCommutator(transpose(L2_is)*L2_is,densityMatrix));
                fmatrix = fmatrix + dissipation1 + dissipation2;
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


