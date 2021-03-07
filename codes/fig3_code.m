%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code generates the plot of the minimum eigenvalue of the Hessian of 
%log-partition function versus inverse temperature beta for the model 
%Hamiltonian H = 10*sum_{i,j} (XiXj + YiYj + ZiZj)+2*sum_i (Xi + Yi + Zi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=5; % size of chain
h=10; % interaction coefficient corresponding to to  XiXj, YiYj, and  ZiZj

betarange=[0:0.001:.10,0.1:0.01:0.85]; % range of inverse temperature
%   to save time the range of beta can be changed to
%   betarange=[0:0.01:.10,0.1:0.01:3];

gap=[]; % minimum eigenvalue is stored in gap
for beta=betarange
    n=6*L; % total number of interaction coefficients
    
    %setting up the coefficients for XiXj, YiYj, and  ZiZj
    coeff=h*ones(3*L,1);
    %setting up the coefficients for Xi, Yi, and  Zi
    coeff(3*L+1:6*L,1)=2*ones(3*L,1);
    
    Iden=eye(n);
    eps=1e-5; % the step size e when evaluating Hessian using (f(x+e)-f(x))/e
    nabla=zeros(n,n); % gradient
    
    % computing the Hessian of log-partition function using differences (f(x+e)-f(x))/e
    for k=1:n
        for j=1:k
            Mean1 = grad(beta,L,coeff); % grad() returns the gradient 
            Mean2 = grad(beta,L,coeff+eps*Iden(:,k)); 
            % computing the Hessian entries:
            nabla(k,j)=(Mean2(j)-Mean1(j))/eps;
            nabla(j,k)=nabla(k,j);
        end
    end
    
    [~,lambda]=eig(nabla); % computing the egienvalues of Hessian
    lambda=sort(diag(lambda));
    gap=[gap real(lambda(1))]; % gap stores the min eigenvalue
end

% fitting a curve on the min eigenvalue vs beta: 
[ff,dc] = fit(betarange', gap', 'a*x^(1.5)*exp(-b*x^2)', 'StartPoint', [1,1]);
plot(ff,betarange,gap)

% generating the figure including the min eigenvalue vs beta and the fitted
% curve
figure;

% Create axes
axes1 = axes('Position',...
    [0.242707117852975 0.355289421157685 0.515248523419927 0.569710578842315]);
hold(axes1,'on');

% Create plot
plot(betarange,gap,'DisplayName','Minimum eigenvalue','Marker','.','LineWidth',2,...
    'Color',[0.635294139385223 0.0784313753247261 0.184313729405403]);
hold on
% Create plot
plot(ff);

% Create ylabel
ylabel('\lambda_{min}(\nabla^2 log Z_{\beta}(\lambda)',...
    'HorizontalAlignment','center',...
    'FontSize',18);

% Create xlabel
xlabel('\beta','FontSize',22);

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 0.3]);
box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontSize',20,'XTick',[0 0.1 0.2 0.3]);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.484325345161964 0.804725823639979 0.272397094430993 0.119251497005988],...
    'FontSize',20,...
    'EdgeColor',[1 1 1]);
