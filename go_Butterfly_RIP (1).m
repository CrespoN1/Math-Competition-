

function go_Butterfly()

%
% Clears any previous plots that are open in MATLAB
clf;


%
% Time Information / Initialization
%
Tstart = 0;            % Simulation starts a tstart (initial value)
Tstop = 30;           % Simulation runs until time = TFinal
         % Simulation starts a tstart (initial value)
           % Simulation runs until time = TFinal



%
% Initial Values
%
Fm0 = 1.125; % Initial Population for Mated Females
Fv0 = 36.38; %initial value for virgin females
M0 = 37.5; % Initial Population for Males
W0= 20.0;%200000;                      % Initial Population for Wasps
E0= 0;
Initial_Values = [Fm0 Fv0 M0 W0 E0]; % Stores initial values in vector


%
% ode45 is matlab's ode solver
%
options=odeset('RelTol',1e-4);
[t,sol] = ode23s(@f,[Tstart Tstop],Initial_Values,options)


%
% storing solutions for each variable, theta_k.
% 
Fm  = sol(:,1);    %gives us Mated Females
Fv =  sol(:,2);    %gives us Virgin Females
M  = sol(:,3);  %gives us Males
W = sol(:,4);   %
E = sol(:,5);


%
% Plotting solutions
%
plot_Time_Evolutions(t,Fm,Fv,M,W,E)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: RHS vector of the problem: this function evaluates the 
%           RHS of the ODEs and passes it back to the ode45 solver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dvdt = f(t,sol)


%
% Components of vectors
%
Fm = sol(1);         % Mated Females
Fv = sol(2);         % Virgin Females
M = sol(3);           % Males
E= sol(4);            % Eggs
W = sol(5);           % Wasps



%
% ODE Parameter Values
%
   


b = 20; %eggs per fm 
dfv = .1; %death rate of virgin females
dm = .125;%.125;%death rate of males and wasps
dw = .15;
dfm = .15; %.15; %enhanced death rate of mated females
de = 0.7; %natural/outside of model death of eggs beed to fix for single day?
omega1 = 0.069; %rate of interaction between mated females and wasps
omega2 = 0.069;
sigma = 0.007;%.0071; %anti aphrodisiac
alpha = .12; %rate of interaction virgin females and males
cfm = 50;
cm = 50;
cw = 34;

%
% ODES (RHS)
%

%With Wasp Population
dEdt = b*dfm*Fm - (omega2 + sigma)*(Fm*W)-de*E;%b*(alpha+sigma)*Fv*M - (omega + sigma)*(Fm*W)-de*E;%b*dfm*Fm - (omega + sigma)*(Fm*W)-de*E;
dFvdt = (E/2)-(alpha+sigma)*(Fv)*(M)-(dfv)*(Fv);
dMdt = (E/2)*(1-(M/cm)) -(dm*M);
dFmdt = (alpha + sigma)*(Fv*M)*(1-(Fm/cfm))-(dfm*Fm);
dWdt = (omega1 + sigma)*(W*Fm)*(1-(W/cw))-dw*W; %(omega + sigma)*(Fm*W)-dmw*W;



%
% Vector to be evaluated
%
dvdt = [dFmdt dFvdt dMdt dWdt dEdt]';



return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots the time evolutions (solutions to ODEs)!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Time_Evolutions(t,Fm,Fv,M,W,E)

%
% Plot Attributes
%
lw = 4;  % LineWidth (how thick the lines should be)
ms = 25; % MarkerSize (how big the plot points should be)
fs = 18; % FontSize (how big the font should be for labels)

%
% PLOT 1: Populations vs. Time
%
figure(1)
plot(t,Fm,'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
plot(t,Fv,'g.-','LineWidth',lw,'MarkerSize',ms); hold on;
plot(t,M,'r.-','LineWidth',lw,'MarkerSize',ms); hold on;
plot(t,W,'m.-','LineWidth',lw,'MarkerSize',ms); hold on;
plot(t,E,'y.-','LineWidth',lw,'MarkerSize',ms); hold on;
hold off;
xlabel('Time')
ylabel('Population sizes (x10,000)')
legend('Mated Females' ,'Virgin Females', 'Males', 'Wasps', 'Eggs')
set(gca,'FontSize',fs)
set(legend,'FontSize',fs)

% PLOT 2: Phase Plane Plot
%
figure(2)
plot(Fm,W,'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Mated f');
ylabel('Wasps');
set(gca,'FontSize',fs);





