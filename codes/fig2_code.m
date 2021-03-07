%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code generates the plot of log-partition function versus interaction
%coefficients μxx and μyy for the model Hamiltonian
%H = sum_{i,j} (μxx XiXj + μyy YiYj + μzz ZiZj)+sum_i (μx Xi + μyYi + μzZi)
%when μzz=2, μx=μy=μz=3.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=0.05;  % beta is the inverse temperature
L=5;    % the number of qubits on the chain
n=6*L;  % the total number of interaction coefficients
h1=2;   % μzz=2
h2=3;   % μx=μy=μz=3

% assigning the coefficients for all terms
coeff=[h1*ones(3*L,1);h2*ones(3*L,1)];

% determining the range of μxx and μyy
h=2;
xx=-5*h:0.2:5*h;
yy=xx;

% creating a mesh containing μxx and μyy
[axx,ayy]=meshgrid(xx,yy);
sizeX=size(xx);
fz=zeros(sizeX(2),sizeX(2));    % fz contains log-partition function values
% obtain the value of log-partition function for all choices of μxx and μyy
for k=1:sizeX(2)
    for l=1:sizeX(2)
        coeff(1:3:3*L,1)=axx(k,l)*ones(size(1:3:3*L,2),1);
        coeff((1:3:3*L)+2,1)=ayy(k,l)*ones(size(1:3:3*L,2),1);
        fz(k,l)=objfunc(beta,L,coeff);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here we generate the 3D plot

% getting the corners of the domain
min_x = min(min(axx));
min_y = min(min(ayy));
max_x = max(max(axx));
max_y = max(max(ayy));

% scaling image between [0, ww] to use a custom color map for it
ww=256;
minplaneimg = min(min(fz)); % find the minimum
scaledimg = (floor(((fz - minplaneimg) ./ ...
    (max(max(fz)) - minplaneimg)) * ww)); % perform scaling


% convert the image to a true color image with the jet colormap.
newmap = brighten(pink(ww),0.2);

colorimg = ind2rgb(scaledimg,newmap);

figure; hold on;

% setting the colormap and generating the plot
newmap2 = brighten(othercolor('OrRd4'),0.);
sc=surf(axx,ayy,fz,'MeshStyle','column','FaceColor','interp','EdgeColor','flat');
hold on
[~,hc] = contourf(axx,ayy,fz);
hc.ContourZLevel = min(min(abs(fz)))*0.95;
hc.LevelStep=(max(max(abs(fz)))-min(min(abs(fz))))*1e-1*.5;
hc.Fill=1;
colormap(flip(othercolor('OrRd7'),1));
view(45,15);

% labels
zlabel('Log Z_{\beta}(\lambda)');
grid('on');
xlabel('J_{xx}')
ylabel('J_{yy}')

