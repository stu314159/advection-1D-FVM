% advec_1d.m

clear
clc
close('all')

x_left = -20;
x_right = 20;


T_end = 10;
N = 1000; % number of cells!!
U = 1;
method = 8;

% Methods

% linear methods
% 1 - upwind
% 2 - Lax-Wendroff
% 3 - Beam-Warming
% 4 - Fromm

% High-resolution limiters
% 5 - minmod
% 6 - superbee
% 7 - MC
% 8 - van Leer

switch (method)
    
    case 1
        phi = @phi_upwind;
        
    case 2
        phi = @lax_wendroff;
        
    case 3
        phi = @beam_warming;
        
    case 4
        phi = @fromm;
        
    case 5
        phi = @minmod_phi;
        
    case 6
        phi = @superbee;
        
    case 7
        phi = @MC_phi;
        
    case 8
        phi = @van_Leer;
        
end



indx = (1:N);
gcoord = zeros(N+1,1);
nodes = zeros(N,2);

% do this just so I have a similar interface for 1D as to 2D.
gcoord(:)=linspace(x_left,x_right,N+1);
nodes(:,1)=indx;
nodes(:,2)=indx+1;
Lk = gcoord(nodes(:,2))-gcoord(nodes(:,1)); dx = gcoord(2)-gcoord(1);
Lc = 0.5*(gcoord(nodes(:,1))+gcoord(nodes(:,2)));
n_ke = ones(N,2); 
n_ke(:,1)=-1; % local mapping is x_l = side 1, x_r = side 2



ElmtToElmnt = zeros(N,2);
ElmntToElmnt(:,1)=circshift(indx,[0 1]);
ElmntToElmnt(:,2)=circshift(indx,[0 -1]);

C = 0.8;
dt = C*min(Lk);

Num_ts = ceil(T_end/dt);
plot_freq = 1;
time_step_rep_freq = 100;

fprintf('Final time = %g, Number of time steps = %d.\n',Num_ts*dt, Num_ts);

% build initial condition
gauss_init = exp(-((Lc + 5).^2)/1);
square_init = zeros(N,1);
sq_ind = find((Lc<7)&(Lc>3));
square_init(sq_ind)=square_init(sq_ind)+1;

q = gauss_init + square_init;
% plot initial condition
plot(Lc,q)

lamda = U*dt/dx;

for ts = 1:Num_ts
    
    if(mod(ts,time_step_rep_freq)==0)
        fprintf('Executing time step number %d.\n',ts);
    end
    
    q = Godunov_1D_general(q,lamda,phi);
    
    if(mod(ts,plot_freq)==0)
        plot(Lc,q)
        grid on
        axis([x_left x_right 0 1.1]);
        drawnow
        
    end
    
    
end

