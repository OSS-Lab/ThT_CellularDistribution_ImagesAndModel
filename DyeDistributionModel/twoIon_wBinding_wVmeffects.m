
%driver script for the ion flows model. 
%This model involves two ions, ThT and TMRM, where flux across plasma and 
%mitochondria membrane is simulated. For ThT, binding in both 
%cytoplasm and mitochondrial conpartments is also simulated where the mitochondrial bound compartment can
%effect the mitochondrial MP. The flux accounts only for passive 
%movement under the influence of MP and conc. 

clear 
clc

%% INITIALISE

%call script that creates all constants and relevant parameters
twoIon_wBinding_wVmeffects_constants

%Set initial ion condition. NOTE: steady state is high outside and zero in plasma. 
y0 = zeros(6,1);

%% RUN SIMS
%Specify the relative and absolute error tolerance:
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

%set simulation time and run 1st simulation
simT = 5000;
int = 0.1
simPoints = 0:int:simT;
[X,Y] =ode23s(@(t,y0)twoIon_wBinding_wVmeffects_ode(t,y0,p,pS,myC),simPoints,y0,options); 

%Record result from previous sim as new initial conditions and introduce perturbation 
y0=Y(end,:);

mP = p(6,1); % save mitochondrial membrane potential for plotting

% Perturbation
pS(8,1) = 1e-2;   
pS(6,1) = 1e-2;

%set simulation time and run 2nd simulation
simT2 = 1080;
simPoints2 = 0:int:simT2;
 
%continue simualation
[X1,Y1] =ode23s(@(t,y0)twoIon_wBinding_wVmeffects_ode(t,y0,p,pS,myC),simPoints2,y0,options); 

%% PLOT RESULTS

%% Plot combined simulation
%X1=X1+X(end);
X2=[X;X1+X(end)];
Y2=[Y;Y1];

figure(1)
subplot(3,1,1)
plot(X2,Y2(:,1),'linewidth',2,'color','black');    
line(X2,Y2(:,3),'linewidth',2,'LineStyle','--','color','k')
xlabel('Time (sec)');
ylabel('Amount [umoles]');
legend({'ThT_c'; 'ThT_c_b_o_u_n_d'});

subplot(3,1,2)
plot(X2,Y2(:,2),'linewidth',2,'color','red'); 
line(X2,Y2(:,4),'linewidth',2,'LineStyle','--','color','red')
xlabel('Time (sec)');
ylabel('Amount [umoles]');
legend({'ThT_m'; 'ThT_m_b_o_u_n_d'});

% Calculate mitochondrial membrane potential 
mP2 = p(6,1); 
e = length(X);
mit_bound_before = Y2(1:e,4);
mitbound_after = Y2(e+1:end,4);
mtVm = (mit_bound_before/mitV) * delta + mP;
mtVm_a =(mitbound_after/mitV) * delta + mP2;
mtVmf = [mtVm;mtVm_a];
 
subplot(3,1,3), plot(X2,mtVmf.*1000,'linewidth',2,'color','blue')
xlabel('Time (sec)');
ylabel('Membrane Potential (mV)');
legend({'mPotmit'});

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18); 
%% Plot normalised results 

sec = 119/int; % last 2 mins of first simulation
last = simT/int; 
start = last - sec; 
Yb = Y(start:last,:);
Yshort = [Yb;Y1];
Xb = X(start:last);
Xshort=[Xb;X1+Xb(end)];

% Normalsie 
first = Yshort(1,:); 
Yn = Yshort./ first; 

% Calculate mitochondrial membrane potential 
e = length(Yshort);
mit_bound_before = Yshort(1:e,4);
mitbound_after = Yshort(e+1:end,4);
mtVm = (mit_bound_before/mitV) * delta + mP;
mtVm_a =(mitbound_after/mitV) * delta + mP2;
mtVmf = [mtVm;mtVm_a];

figure(2)
colororder({'k','k'})
yyaxis left 
plot(Xshort,log(Yn(:,1)),'linewidth',2,'color','b'); 
line(Xshort,log(Yn(:,3)),'linewidth',2,'LineStyle','--','color','b')
line(Xshort,log(Yn(:,2)),'linewidth',2,'LineStyle',':','color','b');  
line(Xshort,log(Yn(:,4)),'linewidth',2,'LineStyle','-.','color','b')
line(Xshort,log(Yn(:,5)),'linewidth',2,'color','r') 
line(Xshort,log(Yn(:,6)),'linewidth',2,'LineStyle','-.','color','r') 
legend({'ThT_c'; 'ThT_c_b_o_u_n_d';'ThT_m'; 'ThT_m_b_o_u_n_d';'TMRM_c'; 'TMRM_m';'ΔΨm'});

xticks(Xshort(1):119:Xshort(end));
xlim([Xshort(1) Xshort(end)]);
xticklabels({'0','2','4','6','8','10','12','14','16','18','20'})
xlabel('Time (mins)');
ylabel('Log Normalised Amount');


% potential normalised 
f = mtVmf(1); 
mtVmfnorm = mtVmf/f;

yyaxis right 
%plot(X2,mtVmf.*1000,'linewidth',2,'color','k')
plot(Xshort,mtVmfnorm,'linewidth',2,'color','k')
xlabel('Time (mins)');
%ylabel('Membrane Potential (mV)');
ylabel('Normalised Membrane Potential');

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18); 
