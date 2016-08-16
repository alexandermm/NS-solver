function vorticityNS
%Name:           Alexander Martinez-Marquese
%Date:           April 19, 2010
%Language used:  M language (MATLAB script)
%Description:       
%Program to solve flow on a tube using vorticity transport, 
%Poisson radial velocity, and Continuity, in cylindrical coordinates.


clear  all;
format long;

%Physical and geometric constants
CLVZ = 5.0;    %m/s
rho  = 8.3e-4; %g/cm^3
R    = 1.0;    %cm
L    = 10.0;   %cm
mu   = 1.4e-4; %g/cm/s
NR   = 21;
NZ   = 201;
numPs  = NR*NZ;
deltat = 1e-3;

%Bi-CGSTAB constants
res     = 1e-8;
maxIter = 100;

[omg,vz,vr,errors] = NS(CLVZ,rho,R,L,mu,NR,NZ,deltat,res,maxIter, 2,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [omg,vz,vr,errors]=NS(CLVZ,rho,R,L,mu,NR,NZ,deltat,res,maxIter,...
                               opt, graphs)
%Solution of incompressible, laminar Navier-Stokes equations in a duct
%Option of convergence critieria given by opt:
%opt = 1: use max % mass difference stopping criteria (default)
%opt = 2 or another number: use % vorticity change stopping criteria


%Calculate constant values 
deltar = R/(NR-1);
deltaz = L/(NZ-1);
numPs  = NR*NZ;
exactOmg = repmat(2*CLVZ*(0:NR-1)'*deltar/R/R,            NZ,1);
exactVz  = repmat(CLVZ*(1.0 - (0:NR-1)'.^2*deltar^2/R^2), NZ,1);

%Solution matrices
if opt == 4
    omg = repmat(2*CLVZ*(0:NR-1)'*deltar/R/R, NZ,1);
else
    omg = zeros(numPs,1);
end
if opt == 3
    vz  = repmat(CLVZ*(1.0 - (0:NR-1)'.^2*deltar^2/R^2), NZ,1);
else
    vz  = zeros(numPs,1);
end
vr   = zeros(numPs,1);
omg2 = zeros(numPs,1);
bVr  = zeros(numPs,1);

%Values and arrays for timesPent
numZs = zeros(NR,1);
k1    = 1:numPs-1;
k2    = 1:numPs-NR;

%Vector used to get rp using index vector
rOrder = repmat(((1:NR)'-1)*deltar, NZ,1);

%Fill A matrix for vr
[As,Aw,Ap,Ae,An] = fillAvr(numPs,NR,NZ, deltar,deltaz,rOrder);

%Actual mass flow
ind    = 1:NR;
v(ind,1) = CLVZ*(1.0 - (ind-1)'.^2*deltar^2/R^2);
ind  = 2:NR;
m1   = 2*pi*rho*deltar*sum(v(ind).*rOrder(ind)) + pi*rho*v(1)*deltar^2/2;

%Iterations
tnum = 0;
criteria = 1;
tol      = 0.01;
omgTol   = 1e-3;
tic
while criteria
    %CALCULATE OMG
    %Interior points
    for z = 2:NZ-1
        ind       = (2:NR-1)' + (z-1)*NR;
        
        %Calculate interior point constants
        cS =  deltat/deltaz*(mu/rho/deltaz + vz(ind)/2);
        cW =  mu*deltat/rho*(1/deltar^2 - 1/2/deltar./rOrder(ind))+...
              deltat/deltar*vr(ind)/2;
        cP = -mu*deltat/rho*(2/deltar^2+2/deltaz^2+1./rOrder(ind).^2)+...
              deltat*vr(ind)./rOrder(ind) + 1;
        cE =  mu*deltat/rho*(1/deltar^2 + 1/2/deltar./rOrder(ind))-...
              deltat/deltar*vr(ind)/2;
        cN =  deltat/deltaz*(mu/rho/deltaz - vz(ind)/2);
        
        omg2(ind) = cS.*omg(ind-NR) + cW.*omg(ind-1) + cP.*omg(ind) + ...
                    cE.*omg(ind+1)  + cN.*omg(ind+NR);
    end
    %Get indices
    r = 2:NR-1;
    z = 1:NZ;
    %West side 
    ind       = 1 + (z-1)*NR;
    omg2(ind) = 0;
    %South side
    ind       = r;
    omg2(ind) = (vr(ind+NR) - vr(ind))/deltaz - ...
                (vz(ind+1)-vz(ind-1))/2/deltar;
    %East side
    ind       =  NR + (z-1)*NR;
    omg2(ind) = -(vz(ind)-vz(ind-1))/deltar;
    %North side
    ind = r + numPs - NR;
    omg2(ind) = omg2(ind-NR);
    %Make omg = omg2
    omg = omg2;
        
    %CALCULATE VR
    %Solve the Poisson equation for radial velocity to get vr
    %Get bVr
    %South and North points
    ind = (1:NR)';
    bVr(ind) = 0;
    ind = (1:NR)' + numPs - NR;
    bVr(ind) = 0;
    %Interior points
    for z = 2:NZ-1
        ind = (2:NR-1)' + (z-1)*NR;
        bVr(ind) = (omg(ind+NR)-omg(ind-NR))/2/deltaz;
    end
    %Solve pentadiagonal system
    %A = diag(As,-NR)+diag(Aw,-1)+diag(Ap)+diag(Ae,1)+diag(An,NR);
    %vr = A\bVr;
    vr = bi_cgstab(As,Aw,Ap,Ae,An,bVr, res,maxIter,numPs,NR,numZs,k1,k2);
    
    %CALCULATE VZ
    %Solve continuity to get vz
    %South side
    ind     = 1:NR;
    vz(ind) = CLVZ*(1.0 - (ind-1).^2*deltar^2/R^2);
    %Interior points
    for z = 2:NZ
        ind = (2:NR-1)' + (z-1)*NR;
        vz(ind) = -deltaz*vr(ind)./rOrder(ind)- ...
                   deltaz*(vr(ind+1)-vr(ind-1))/2/deltar + vz(ind-NR);
    end
    %North side
    ind = (2:NR-1)' + numPs - NR;
    vz(ind) = vz(ind-NR);
    %West side (east side already zero)
    z = 2:NZ;
    ind = 1 + (z-1)*NR;
    vz(ind) = vz(ind+1);
    %vz
    %Increase time step
    tnum = tnum+1;
    
    if opt == 1
        %Use % mass flow error if opt = 1
        %Calculate mass flow rates and % error
        for z = 1:NZ;
            ind  = (2:NR-1)' + (z-1)*NR;
            ind1 = 1 + (z-1)*NR;
            m(z) = pi*rho*(2*deltar*sum(vz(ind).*rOrder(ind))+vr(ind1)*deltar^2/2);
        end
        %Calculate errors
        perError        = max(abs(m-m1))/m1;
        errors(tnum)    = perError;
        avgErrors(tnum) = mean(abs(m-m1))/m1;
        %Evaluate error criteria
        criteria = (perError > tol) || (tnum < 10);
    end
    if opt == 2
        %Otherwise use % vorticity change error
        %Calculate mass flow rates and % error
        for z = 1:NZ;
            ind  = (2:NR-1)' + (z-1)*NR;
            ind1 = 1 + (z-1)*NR;
            m(z) = pi*rho*(2*deltar*sum(vz(ind).*rOrder(ind))+vr(ind1)*deltar^2/2);
        end
        %Calculate errors
        perError        = max(abs(m-m1))/m1;
        errors(tnum)    = perError;
        avgErrors(tnum) = mean(abs(m-m1))/m1;
        
        if tnum == 1
            perError    = 1;
        else
            perError    = max(abs(omg-oldOmg))/max(abs(omg));
        end
        oldOmg          = omg;
        errors2(tnum)   = perError;
        %Evaluate error criteria
        criteria = (perError > omgTol) || (tnum < 10);
    end
    if opt == 3
        perError = max(abs(exactOmg-omg))/max(abs(exactOmg));
        %Declare errors to send something out of function
        errors = 0;
        %Evaluate error criteria
        criteria = (perError > 0.12882 ) || (tnum < 10);
    end
    if opt == 4
        perError = max(abs(exactVz -vz ))/max(abs(exactVz ));
        %Declare errors to send something out of function
        errors = 0;
        %Evaluate error criteria
        criteria = (perError > 0.098997) || (tnum < 10);
    end
end
elapsedTime = toc;
%Display processing time
disp(sprintf('\n  Processing time: %f', elapsedTime))

%Change error vector to NR*NZ x 1 
errors = errors';

%Display the number of iterations
disp(sprintf('\n  Number of iterations required for convergence: %d',tnum))

if (opt == 1) && (graphs == 1)
    %Plot opt 1 errors
    figure
    plot(errors)
    xlabel('Iteration number','FontSize',12)
    ylabel('Max % mass difference','FontSize',12)
    title('Max % mass difference vs. iteration number','FontSize',12)
    figure
    plot(avgErrors)
    xlabel('Iteration number','FontSize',12)
    ylabel('Avg % mass difference','FontSize',12)
    title('Avg % mass difference vs. iteration number','FontSize',12)
end
if (opt == 2) && (graphs == 1)
    %Plot opt 2 errors (%mass difference and %vorticity change)
    figure
    plot(errors)
    xlabel('Iteration number','FontSize',12)
    ylabel('Max % mass difference','FontSize',12)
    title('Max % mass difference vs. iteration number','FontSize',12)
    figure
    plot(avgErrors)
    xlabel('Iteration number','FontSize',12)
    ylabel('Avg % mass difference','FontSize',12)
    title('Avg % mass difference vs. iteration number','FontSize',12)
    
    figure
    semilogy(errors2)
    xlabel('Iteration number','FontSize',12)
    ylabel('Log of max % omg difference','FontSize',12)
    title('Max % omg difference vs. iteration number','FontSize',12)
end

if graphs == 1
    %Display velocities
    zOrder = [];
    for z = 0:NZ-1
        zOrder = [zOrder; repmat(z*deltaz, NR,1)];
    end
    %Graph portion of points
    ind = [];
    for z = 1:4:NZ
        ind = [ind; (1:4:NR)' + (z-1)*NR];
    end
    z  = zOrder(ind);
    r  = rOrder(ind);
    vZ = vz(ind);
    vR = vr(ind);
    figure
    quivers(z,r,vZ,vR,0.4,1,'cm/s','b')
    xlabel('Vertical distance from bottom of pipe section (cm)','FontSize',12)
    ylabel('Dist. from center of pipe (cm)','FontSize',12)
    title('Velocity plot in pipe after convergence','FontSize',12)
    %Display vorticity
    r = ((1:NR)-1)'*deltar;
    z = ((1:NZ)-1)'*deltaz;
    realOmg     = reshape(omg,NR,NZ);
    figure
    [c,h] = contourf(z,r,realOmg);
    colormap jet;
    cbar_handle = colorbar('location','eastoutside');
    set(get(cbar_handle,'ylabel'),'string','Vorticity (1/s)','FontSize',10)
    daspect([1 1 1])
    set(h,'EdgeColor','none') 
    xlabel('Vertical distance from bottom of pipe section (cm)','FontSize',12)
    ylabel('Dist. from center of pipe (cm)','FontSize',12)
    title('Vorticity plot in pipe after convergence','FontSize',12)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s,w,p,e,n] = fillAvr(numPs,NR,NZ, deltar,deltaz,rOrder)


%Declarations
s = zeros(numPs-NR, 1);
w = zeros(numPs-1 , 1);
p = zeros(numPs   , 1);
e = zeros(numPs-1 , 1);
n = zeros(numPs-NR, 1);

%Fill constant vectors
%Get indices
r = 2:NR-1;
z = 1:NZ;
%West side 
ind = 1 + (z-1)*NR;
p(ind   ) = 1;
%South side
ind = r;
p(ind  ) = 1;
%East side
ind = NR + (z-1)*NR;
p(ind   ) = 1;
%North side
ind = r + numPs - NR;
s(ind-NR) = -1;
p(ind   ) =  1;

%Interior points
for z = 2:NZ-1
    ind = (2:NR-1)' + (z-1)*NR;
    s(ind-NR) =  1/deltaz^2;
    w(ind-1 ) =  1/deltar^2 - 1/2/deltar./rOrder(ind);
    p(ind   ) = -2/deltar^2 - 2/deltaz^2 - 1./rOrder(ind).^2;
    e(ind   ) =  1/deltar^2 + 1/2/deltar./rOrder(ind);
    n(ind   ) =  1/deltaz^2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = bi_cgstab(d, e,f,g, h,  b,  res,maxiter,N,sqN,numZs,k1,k2)
%Description:       
%Bi-STAB (BiConjugate Gradient Stabilized method with preconditionng) 
%algorithm from the paper:
%BI-CGSTAB: A fast and smoothly converging variant of 
%BI-CG for the solution of nonsymmetric linear systems
%by H.A. Van Der Vorst
%Using x0 = b./d and conditioning matrix K = diag(A) = f
%Made to handle pentadiagonal systems using vectorization


%Initial guess
x  = b./f; 
%Get initial error r
r  = b - timesPent(d,e,f,g,h,x,N,sqN,numZs,k1,k2); 
r0 = r;

%Initialize variables
rho = 1; alpha = 1; omega = 1;
v = zeros(N,1); p = zeros(N,1);
iter = 0;

%Actual iterations
while (norm(b-timesPent(d,e,f,g,h,x,N,sqN,numZs,k1,k2),2)/N > res) ...
       && (iter < maxiter)
    lastRho = rho;
    rho     = r0.'*r;
    beta  = (rho/lastRho)*(alpha/omega);
    p     = r + beta*(p - omega*v);

    y     = p./f; 
    v     = timesPent(d,e,f,g,h,y,N,sqN,numZs,k1,k2);
    alpha = rho/(r0.'*v);
    s     = r - alpha*v;
    
    z     = s./f; 
    t     = timesPent(d,e,f,g,h,z,N,sqN,numZs,k1,k2); 
    Kinv  = 1./f; 
    omega = ((Kinv.*t).'*(Kinv.*s))/((Kinv.*t).'*(Kinv.*t));  
    x     = x + alpha*y + omega*z;
    
    r = s - omega*t;
    iter = iter + 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = timesPent(d, e,f,g, h,  x,N,sqN,numZs,k1,k2)
%Description:       
%Function for fast multiplication of pentadiagonal matrices and vectors

b = f.*x + [0;e.*x(k1)] + [g.*x(k1+1);0] + [numZs; d.*x(k2)] ...
    + [h.*x(k2+sqN); numZs];