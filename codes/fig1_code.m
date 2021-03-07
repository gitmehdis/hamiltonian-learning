%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code implements the Gradient Descent method to solve the optimizaiton
%program corresponding to the dual of the maximum entropy optimizaiton. 
%The output is the plot of error in the interaction coefficients versus
%inverse temperature for a generic 2-local Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normdifff=[]; % ell_2 error of the optimum compared to the real coefficients
errparam=1e-3; % standard deviation of the statistical errors 
h=2; % norm of the interaction coefficients
L=5; % size of the chain
n=6*L; % total number of interaction coefficients
coeff=h*(2*rand(n,1)-1); % generating all the interaction coefficients
%coeff=coeff2'; % uncomment to include the parameters used in the paper
% coeff(1:3*L,1)=3*ones(3*L,1);
% coeff(3*L+1:6*L,1)=2*ones(3*L,1);

% range of inverse temperature. To test the code this can be limited to
% shorten the running time
betarange=[[0.01:0.01:0.1 0.15:0.05:0.6] 0.65:0.05:0.85];


for ss=1:10 % number of times local expectations are drawn at random
    normdiff=[]; % the error in ell_2 norm stored here
    
    for beta=betarange
        
        x=coeff;% the initial starting point
        
        % Armijo stepsize rule parameters
        sigma = 0.1;
        bb = 0.5;
         % the local expectations drawn at random:
        lerr=lambdaerr(beta,L,coeff,errparam); 
        % objective function of the dual program is calculated:
        obj=objfunc(beta,L,x)+beta*sum(x.*lerr);
        % gradient of the objective function is calculated:
        g=grad(beta,L,x)+beta*lerr;
        k=0;                    % k = # iterations
        nf=1;					% nf = # function eval.
        
        % gradient descent starts here
        while  (norm(g) > 1e-5 && k<7000) %continue as long as the gradient is large
            d = -g; % steepest descent direction
            a = 1;
            newx = x + a*d;         
            newobj = objfunc(beta,L,newx)+beta*sum(newx.*lerr); % new objective function point
            nf = nf+1;
            rf=0;
            % the size of each step is set here:
            while ( (newobj-obj)/a > sigma*g'*d & rf<10^2) 
                ans1=sigma*g'*d;
                ans2=(newobj-obj)/a;
                a = a*bb;
                newx= x + a*d;
                newobj = objfunc(beta,L,newx)+beta*sum(newx.*lerr);
                nf = nf+1;
                rf=rf+1;
            end
            % progress is printed
            if (mod(k,100)==1) 
                fprintf('%5.0f %5.0f %12.5e %12.5e \n',k,nf,norm(x-coeff)/max(coeff),norm(g));
            end
            x = x + a*d;
            obj=newobj;
            g=grad(beta,L,x)+beta*lerr; %step in the direction of the gradient
            k = k + 1;
        end
        normdiff=[normdiff norm(x-coeff)]; % ell_2 error in the coefficients
    end
    % ell_2 errors for various beta for one random draw is stored
    normdifff=[normdifff normdiff']; 
end
% mean and variance of at each beta is obtained:
meandiff=mean(normdifff')';
vardiff=var(normdifff')';
stdiff=std(normdifff')';

% plotting the data
figure
fill([betarange';flipud(betarange')],[meandiff-stdiff;flipud(meandiff+stdiff)],[0.937254905700684 0.866666674613953 0.866666674613953],'linestyle','none');
line(betarange',meandiff,'LineWidth',3,...
    'Color',[0.600000023841858 0.200000002980232 0])
xlim([0 0.85]);
 

