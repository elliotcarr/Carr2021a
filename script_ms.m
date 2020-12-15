close all, 
clc, clear all

% save_figs = true;
save_figs = false;
addpath('./Carr2020a-master/') % Download from https://github.com/elliotcarr/Carr2020a

font_size = 28;
if save_figs
    path_name = '../../Paper/Figures/';
    addpath('../export_fig-master')
end

Problem = 'A'; breakthrough = false;
% Problem = 'B'; breakthrough = false;
% Problem = 'C'; breakthrough = false;
% Problem = 'D'; breakthrough = false;
% Problem = 'D'; breakthrough = true;
[z,w,Np] =  poles_residues;
N = 10001; % number of nodes for finite volume method (benchmark solution)

%% Test Problems
if strcmp(Problem,'A')
    n = 4; % number of species
    L = 100; % length of medium
    m = 2; % number of layers
    l = [50]; plot_layers = []; % interfaces
    D = [0.3,0.3]; % D(i) - diffusivity in layer i
    v = [0.2,0.2]; % v(i) - velocity in layer i
    theta = [1.0,1.0]; % theta(i) - volumetric water content in layer i  
    Rflag = 'layers'; Rlayers = ones(m,1); % R(i) - retardation factor in layer i
    f = zeros(m,n); % f(i,j) - initial concentration in layer i for species j
    gamma = zeros(m,n); % gamma(i,j) - zero order production in layer i for species j
    xticks = 0:20:100; % tick marks along x-axis for concentration plots
    
    % Reaction matrices
    % M(i,k) - rate constant for first-order production/decay of species i OR first-order production
    % of speces i from species k.
    M = [-0.05 0 0 0; 0.05 -0.02 0 0; 0 0.02 -0.01 0; 0 0 0.01 -0.005];
    
    % Time
    tspan = [200,400];
    
    % Inlet boundary condition
    a0 = 1.0; b0 = 0.0;
    inlet{1} = {a0, b0, @(t) 1, @(s) 1/s}; % species 1
    inlet{2} = {a0, b0, @(t) 0, @(s) 0}; % species 2
    inlet{3} = {a0, b0, @(t) 0, @(s) 0}; % species 3
    inlet{4} = {a0, b0, @(t) 0, @(s) 0}; % species 4
  
    % Outlet boundary condition
    aL = 0; bL = 1; gL = @(t) 0; GL = @(s) 0;
    outlet{1} = {aL, bL, gL, GL};
    outlet{2} = {aL, bL, gL, GL};
    outlet{3} = {aL, bL, gL, GL};
    outlet{4} = {aL, bL, gL, GL};
    
elseif strcmp(Problem,'B')
    n = 4; % number of species
    L = 1; % length of medium
    m = 3; % number of layers
    l = [0.3,0.5]; plot_layers = l; % interfaces
    D = [0.01,0.004,0.01]; % D(i) - diffusivity in layer i
    v = [1.0,0.24,1.0]; % v(i) - velocity in layer i
    theta = [0.12,0.5,0.12]; % theta(i) - volumetric water content in layer i 
    Rflag = 'layers'; Rlayers = ones(m,1); % R(i) - retardation factor in layer i
    f = zeros(m,n); % f(i,j) - initial concentration in layer i for species j
    gamma = zeros(m,n); % gamma(i,j) - zero order production in layer i for species j
    xticks = [0:0.2:1]; % tick marks along x-axis for concentration plots
    
    % Reaction matrices
    % M(i,k) - rate constant for first-order production/decay of species i OR first-order production
    % of speces i from species k.
    M = [-0.75 0 0 0; 0.75 -0.5 0 0; 0 0.5 -0.2 0; 0 0 0.2 -0.1];

    % Time
    tspan = [0.6,1];
    
    % Inlet boundary condition
    a0 = theta(1)*v(1); b0 = theta(1)*D(1); g0 = @(t) theta(1)*v(1)*1.0; G0 = @(s) theta(1)*v(1)*1.0/s;
    inlet{1} = {a0, b0, g0, G0}; % species 1
    inlet{2} = {a0, b0, @(t) 0, @(s) 0}; % species 2
    inlet{3} = {a0, b0, @(t) 0, @(s) 0}; % species 3
    inlet{4} = {a0, b0, @(t) 0, @(s) 0}; % species 4
  
    % Outlet boundary condition
    aL = 0; bL = 1; gL = @(t) 0; GL = @(s) 0;
    outlet{1} = {aL, bL, gL, GL};
    outlet{2} = {aL, bL, gL, GL};
    outlet{3} = {aL, bL, gL, GL};
    outlet{4} = {aL, bL, gL, GL};

elseif strcmp(Problem,'C')
    n = 4; % number of species
    L = 40; % length of medium
    m = 5; % number of layers
    l = [10,15,20,22]; plot_layers = l; % interfaces
    D = [0.08,0.008,0.08,0.008,0.08]; % D(i) - diffusivity in layer i
    v = [0.4,0.04,0.4,0.04,0.4]; % v(i) - velocity in layer i
    theta = [0.09,0.9,0.09,0.9,0.09]; % theta(i) - volumetric water content in layer i
    Rflag = 'layers'; Rlayers = [1.0,0.8,1.0,0.8,1.0]; % R(i) - retardation factor in layer i
    f = zeros(m,n); % f(i,j) - initial concentration in layer i for species j
    gamma = zeros(m,n); % gamma(i,j) - zero order production in layer i for species j
    xticks = [0:10:40]; % tick marks along x-axis for concentration plots
    
    % Reaction matrices
    % M(i,k) - rate constant for first-order production/decay of species i OR first-order production
    % of speces i from species k.
    M = [-0.075,0,0,0; 0.75*0.075 -0.05 0.2*0.02 0; 0.25*0.075 0.5*0.05 -0.02 0; 0 0.5*0.05 0.8*0.02 -0.045];   

    % Time
    tspan = [100,400];
    
    % Inlet boundary condition
    t0 = 200; beta = -0.01;
    a0 = 1.0; b0 = 0.0; g0 = @(t) 1.0-exp(beta*t); G0 = @(s) 1.0/s-1.0/(s-beta);  
    inlet{1} = {a0, b0, g0, G0}; % species 1
    inlet{2} = {a0, b0, @(t) 0, @(s) 0}; % species 2
    inlet{3} = {a0, b0, @(t) 0, @(s) 0}; % species 3
    inlet{4} = {a0, b0, @(t) 0, @(s) 0}; % species 4  
  
    % Outlet boundary condition
    aL = 0; bL = 1; gL = @(t) 0; GL = @(s) 0;
    outlet{1} = {aL, bL, gL, GL};
    outlet{2} = {aL, bL, gL, GL};
    outlet{3} = {aL, bL, gL, GL};
    outlet{4} = {aL, bL, gL, GL};
    
elseif strcmp(Problem,'D')
    n = 4; % number of species
    L = 40; % length of medium
    m = 5; % number of layers
    l = [10,15,20,22]; plot_layers = l; % interfaces
    D = [0.08,0.008,0.08,0.008,0.08]; % D(i) - diffusivity in layer i
    v = [0.4,0.04,0.4,0.04,0.4]; % v(i) - velocity in layer i
    theta = [0.09,0.9,0.09,0.9,0.09]; % theta(i) - volumetric water content in layer i 
    Rflag = 'species'; Rspecies = [1.0,0.9,0.5,0.8]; % R(j) - retardation factor for species j
    f = zeros(m,n); % f(i,j) - initial concentration in layer i for species j
    gamma = zeros(m,n); gamma(4,1) = 0.01; % gamma(i,j) - zero order production in layer i for species j
    xticks = [0:10:40]; % tick marks along x-axis for concentration plots

    % Reaction matrix
    % M(i,k) - rate constant for first-order production/decay of species i OR first-order production
    % of speces i from species k.
    M = [-0.075,0,0,0; 0.6*0.075 -0.05 0.2*0.02 1.0*0.045; 0.25*0.075 0.5*0.05 -0.02 0; 0.15*0.075 0.5*0.05 0.8*0.02 -0.045];

    % Time
    if breakthrough
        tspan = linspace(0,400,41); tspan(1) = 1e-2;
    else
        tspan = [100,400];
    end
    
    % Inlet boundary condition
    a = pi/(800);
    a0 = 1.0; b0 = 0.0; g0 = @(t) 0.5*(1+cos(a*t)); G0 = @(s) 0.5*(1/s + s/(s^2+a^2));
    inlet{1} = {a0, b0, g0, G0}; % species 1
    inlet{2} = {a0, b0, @(t) 0, @(s) 0}; % species 2
    inlet{3} = {a0, b0, @(t) 0, @(s) 0}; % species 3
    inlet{4} = {a0, b0, @(t) 0, @(s) 0}; % species 4
  
    % Outlet boundary condition
    aL = 0; bL = 1; gL = @(t) 0; GL = @(s) 0;
    outlet{1} = {aL, bL, gL, GL};
    outlet{2} = {aL, bL, gL, GL};
    outlet{3} = {aL, bL, gL, GL};
    outlet{4} = {aL, bL, gL, GL};
    
end
l = [l,L];
x = linspace(0,L,N); % fvm nodes
% semi-analytical/plotting nodes
if breakthrough
    xp = L;
else
    xp = linspace(0,L,501);
end

%% Semi-analytical solution
tic
% Case I: Non-equal retardation factors across layers
if strcmp(Rflag,'layers')
    
    R = repmat(reshape(Rlayers,m,1),1,n);
    
    % Apply species-decoupling transformation
    [Q,Lambda] = eig(M);
    lambda = diag(Lambda);
    inlett = cell(n,1);
    outlett = cell(n,1);
%     ut = Q\eye(n,1);
%     for j = 1:n
%         inlett{j} = {a0, b0, @(t) ut(j)*inlet{1}{3}(t), @(s) ut(j)*inlet{1}{4}(s)};
%     end
    Qinv = Q\eye(n);
    for j = 1:n
        inlett{j} = cell(4,1);
        inlett{j}{1} = a0;
        inlett{j}{2} = b0;
        inlett{j}{3} = @(t) g0t_func(j,t,inlet,Qinv,n);
        inlett{j}{4} = @(s) G0s_func(j,s,inlet,Qinv,n);
        outlett{j} = cell(4,1);
        outlett{j}{1} = aL;
        outlett{j}{2} = bL;
        outlett{j}{3} = @(t) g0t_func(j,t,outlet,Qinv,n);
        outlett{j}{4} = @(s) G0s_func(j,s,outlet,Qinv,n);        
    end
    ft = Q\f';
    gammat = Q\gamma';    
    
    % Semi-analytical solution using Laplace transform
    ct = zeros(length(xp),length(tspan),n);
    for j = 1:n % Loop over number of species
        for p = 1:length(tspan)
            t = tspan(p);
            for k = 1:length(xp)
                Ct = @(s) Cfunc(s,xp(k),v,D,theta,-lambda(j)*ones(1,m),Rlayers,gammat(j,:),l,ft(j,:),inlett{j},outlett{j});
                ct(k,p,j) = inverse_laplace_transform(Ct,t,Np,w,z);
            end
        end
    end
    
    % Reverse species-decoupling transformation
    c = zeros(size(ct));
    for i = 1:n
        for j = 1:n
            c(:,:,i) = c(:,:,i) + Q(i,j)*ct(:,:,j);
        end
    end

% Case II: Non-equal retardation factors across species
elseif strcmp(Rflag,'species')
    
    R = repmat(reshape(Rspecies,1,n),m,1);
    
    % Semi-analytical solution using Laplace transform
    c = zeros(length(xp),length(tspan),n);
    for p = 1:length(tspan)
        t = tspan(p);
        for k = 1:length(xp)
            C = @(s) Cfunc_wrap(s,xp(k),v,D,theta,M,Rspecies,gamma,l,f,inlet,outlet,n,m);
            % Inverse Laplace transform
            invL = zeros(n,1);
            for r = 1:2:Np
                s = z(r)/t;
                invL = invL-w(r)*C(s);
            end
            c(k,p,:) = 2*real(invL)/t;
        end
    end
end
toc

%% Numerical solution
tic
cnt = numerical_solution_ms(tspan,R,D,v,M,gamma,theta,l,f,N,m,inlet,outlet,n);
toc
cnt = cnt'; cn = zeros(length(x),length(tspan),n);
for j = 1:n
    if length(tspan) == 1
        cntp = cnt(:,end);
        cn(:,1,j) = cntp((j-1)*N+1:j*N);
    else
        for p = 1:length(tspan)
            cntp = cnt(:,p);
            cn(:,p,j) = cntp((j-1)*N+1:j*N);
        end
    end
end

%% Compute maximum absolute errors (Table 3 in manuscript)
ix = (N-1)/(501-1);
err = max(abs(c-cn(1:ix:N,:,:)),[],'all');
fprintf('\n')
fprintf(['Entry in Table 3 for Problem ',Problem,'\n'])
fprintf('\\num{%1.1e}\n',err);
xscal = 1; tscal = 1;

%% Plot solutions (Figure 2 in manuscript)
markers = {'o','s','^','v'};
colors = [127,205,187; 65,182,196; 29,145,192; 34,94,168]/255;

if ~breakthrough
    for p = 1:length(tspan)
        figure;
        vec = get(gcf,'Position');
        set(gcf,'Color','w','Position',vec.*[1,1,1.3,1])
        leg_str = cell(n,1);
        h = zeros(n,1);
        for j = 1:n            
            if ~strcmp(Problem,'A')
                for i = 1:length(plot_layers)
                    plot([plot_layers(i),plot_layers(i)],[0,1],'-','Color',0.7*ones(3,1))
                    hold on
                end
            end
            box on            
            plot(xp*xscal,c(:,p,j),'Color',colors(j,:),'Linewidth',2.0);
            hold on
            plot(x*xscal,cn(:,p,j),['-',markers{j}],'Color',colors(j,:),'MarkerSize',12,...
                'MarkerIndices',round(linspace(1,N,41)),'MarkerFaceColor',colors(j,:));
            h(j) = plot(x*xscal,NaN*cn(:,p,j),['-',markers{j}],'Color',colors(j,:),'MarkerSize',12,...
                'MarkerIndices',round(linspace(1,N,41)),'MarkerFaceColor',colors(j,:),'Linewidth',2.0); % dummy for legend
            axis([0 L 0 1])
            set(gca,'Fontsize',font_size-2,'TickLabelInterpreter','latex','Xtick',xticks,'XMinorTick','on',...
                'YTick',[0,0.5,1],'YMinorTick','on','TickDir','out','Color',0.95*ones(1,3))
            xlabel('$x$ [$\textrm{m}$]','Interpreter','LaTeX','Fontsize',font_size)
            ylabel('$c_{j}(x,t)$ [$\textrm{mol}\,\textrm{L}^{-1}$]','Interpreter','LaTeX','Fontsize',font_size)
            text(0.6*L,0.8,['Problem ',Problem],'Interpreter','LaTeX','Fontsize',font_size-4)
            if tspan(p)*tscal <= 1
                text(0.6*L,0.72,['$t = ',num2str(tspan(p)*tscal),'\,\textrm{day}$'],'Interpreter','LaTeX','Fontsize',font_size-4)
            else
                text(0.6*L,0.72,['$t = ',num2str(tspan(p)*tscal),'\,\textrm{days}$'],'Interpreter','LaTeX','Fontsize',font_size-4)
            end
            leg_str{j} = ['$',num2str(j),'$'];
        end
        leg = legend(h,leg_str,'Interpreter','LaTeX','Fontsize',font_size-8,'Location','NorthEast');
        leg.Title.String = 'Species';
        drawnow
        if save_figs
            feval('export_fig',[path_name,'Problem',Problem,'_Time',num2str(p)],'-pdf')
        end
    end
end

%% Plot breakthrough curves for Problem D (Figure 3 in manuscript)
if breakthrough
    figure
    vec = get(gcf,'Position');
    set(gcf,'Color','w','Position',vec.*[1,1,1.3,1])
    for j = 1:n
        plot(tspan,c(end,:,j),'Color',colors(j,:),'Linewidth',2.0)
        hold on
        plot(tspan,cn(end,:,j),['-',markers{j}],'Color',colors(j,:),'MarkerSize',12,...
            'MarkerFaceColor',colors(j,:));
        h(j) = plot(tspan,NaN*cn(end,:,j),['-',markers{j}],'Color',colors(j,:),'MarkerSize',12,...
            'MarkerFaceColor',colors(j,:),'Linewidth',2.0);
        axis([0 400 0 0.6])
        set(gca,'Fontsize',font_size-6,'TickLabelInterpreter','latex','XMinorTick','on',...
            'YTick',[0:0.1:0.6],'YMinorTick','on','TickDir','out','Color',0.95*ones(1,3))
        xlabel('$t$ $[\textrm{days}]$','Interpreter','LaTeX','Fontsize',font_size-4)
        ylabel('$c_{j}(L,t)$ [$\textrm{mol}\,\textrm{L}^{-1}$]','Interpreter','LaTeX','Fontsize',font_size-4)
        leg_str{j} = ['$',num2str(j),'$'];
    end
    leg = legend(h,leg_str,'Interpreter','LaTeX','Fontsize',font_size-8,'Location','NorthWest');
    leg.Title.String = 'Species';
    drawnow
    if save_figs
        feval('export_fig',[path_name,'Problem',Problem,'_Breakthrough'],'-pdf')
    end
end