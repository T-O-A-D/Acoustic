% This program describes a moving 1-D wave
% using the 2-D finite difference method
clc
close all;
clear all;
%-------------------------------------------------------------------------%
%% Initialization
Nx = 101;       % Number of grids/steps in x 
dx = 1;         % Step size of x
x(:,1) = (0:Nx-1)*dx;   % x-axis
          
T = 501;        % Total number of time steps
f = 10;         % Frequency of source
dt = 0.001;     % size of Time-Step
t(:,1)= (0:T-1)*dt;     % Time-axis
v = 500;        % Wave velocity
c = v*(dt/dx);   % CFL condition, Try to keep < 1
U = zeros(T,Nx);  % U(x,t) = U(time,space)
s1 = floor(T/f);  
%-------------------------------------------------------------------------%
%% Initial condition
U((1:s1),1) = sin(2*pi*f.*t(1:s1));
U((1:s1),2) = sin(2*pi*f.*t(1:s1));
V = U;      % V is another wave with identical initial condition
W = U;
%-------------------------------------------------------------------------%
%% Finite Difference Scheme
for j = 3:T-1
    for i = 2:Nx-1
        
        %% Clamped End
        U1 = 2*U(j-1,i)-U(j-2,i); %finite difference in time
        U2 = U(j-1,i-1) - 2*U(j-1,i) + U(j-1,i+1); %finite difference in space
        U(j,i) = U1 + c*c.*U2; 
        
        %% Simply Support
        V1 = 2*V(j-1,i)-V(j-2,i);%finite difference in time
        V2 = V(j-1,i-1)-2*V(j-1,i)+V(j-1,i+1); %finite difference in space
        V(j,i) = V1 + c*c.*V2;

        %% Free End
        W1 = 2*W(j-1,i)-W(j-2,i); %finite difference in time
        W2 = W(j-1,i-1)-2*W(j-1,i)+W(j-1,i+1); %finite difference in space
        W(j,i) = W1 + c*c.*W2;
    end     
    
     % Artificial & approximate Non-reflecting Boundary 
     % no end
     % U(j+1,Nx) = 0.5*( U(j,Nx)+U(j,Nx-1)); 

     W(j+1,Nx) = W1; % Free End
     V(j+1,Nx) = 0; % Simply Support
     U(j+1,Nx) = 0; % Clamped end
     U(j,Nx-1) = 0; % Clamped end
end
%-------------------------------------------------------------------------%
%% Plot for the travelling waves
for j = 1:T   
  
  %% case 1 - Clamped End 
  subplot 311  
  plot(x,U(j,:),'r','linewidth',2);
  line([Nx Nx],[-1 1],'Color','k','LineStyle','--','LineWidth',3);
  text(70,0.8,'Clamped End --->','FontSize',14)
  grid on;
  axis([min(x) max(x)+3 -1.25 1.25]);
  %xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);
  title1 = ['TIME STEP = ',num2str(j), '   TIME = ',num2str(t(j)),' second'];
  title(title1 ,'fontsize',14);
  h=gca; 
  get(h,'FontSize'); 
  set(h,'FontSize',14);
  
  %% case 2 - Free End 
  subplot 312
  plot(x,W(j,:),'g', 'linewidth',2);
  line([Nx Nx],[-1 1],'Color','k','LineStyle','--','LineWidth',3);
  text(70,0.8,'Free End --->','FontSize',14)
  grid on;
  axis([min(x) max(x)+3 -1.25 1.25]);
  %xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14); 
  %title3 = ['TIME STEP = ',num2str(j), '   TIME = ',num2str(t(j)),' second'];
  %title(title3 ,'fontsize',14);                            
  h = gca; 
  get(h,'FontSize'); 
  set(h,'FontSize',14);
  
  %% case 3 - Simply Supported End
  subplot 313 
  plot(x,V(j,:),'b','linewidth',2);
  line([Nx Nx],[-1 1],'Color','k','LineStyle','--','LineWidth',3);
  text(70,0.8,'Simply Supported --->','FontSize',14)
  grid on;
  axis([min(x) max(x)+3 -1.25 1.25]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14); 
  %title2 = ['TIME STEP = ',num2str(j), '   TIME = ',num2str(t(j)),' second'];
  %title(title2 ,'fontsize',14);                            
  h = gca; 
  get(h,'FontSize'); 
  set(h,'FontSize',14);
  
  fh = figure(1);
  set(fh, 'color', 'white'); 
end
%-------------------------------------------------------------------------%
