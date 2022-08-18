
% Thomas Pile
% MSc Control and Robotics, Sheffield Hallam
% LMI code from Jianchi Chen, Uni of Leicester (~2010)
% LMI solution for coprime factorisation from Prempain, also UoL
% 23 May 2018

% This code uses LMIs and an algorithm by Prempain to produce a left
% coprime factorisation of a plant by solving an H2 optimisation problem.

%%
% 1. Make the GS model using LPV LMI model

% This is from my proof of concept controller 1
A = ...
        [0         0         0    0.9900         0         0;
         0         0         0         0    0.9900         0;
         0         0         0         0         0    0.9900;
         0         0         0   -7.1429    0.6071    0.6071;
         0         0         0   -0.6071   -7.1429   -0.6071;
         0         0         0         0         0  -18.1818];
B = ...
        [0         0         0;         
         0         0         0;         
         0         0         0;       
  109.2857         0         0;        
         0  109.2857         0;         
         0         0  278.1818];
C = ...
  [57.2958         0         0         0         0         0;
         0   57.2958         0         0         0         0;
         0         0   57.2958         0         0         0];
D = ...
    [0     0     0;
     0     0     0;
     0     0     0];
Zmin=1; Zmax=4;
pv = pvec('box',[Zmin Zmax]);
% affine model:
s0 = ltisys(A,B,C,D);
s1 = ltisys([0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0.9 0 0; 
             0 0 0 0 0.7 0; 
             0 0 0 0 0 0.5],...
             zeros(6,3),zeros(3,6),zeros(3,3),0); 
pdG = psys(pv,[s0 s1]);
tpG = aff2pol(pdG);
minfo(tpG)

%%
% 2. Extract the linear models at grid points

% sizes of the main matrices
szA = size(A);
szB = size(B); % 6 3   
szC = size(C); % 3 6  
% This needs to be updated to use the counting method later
numverts = 2;
A = zeros(szA(1),szA(2),numverts);
% extract B, C, D
pdsi = psinfo(tpG,'sys',1);
B = pdsi(1:szB(1), szA(2)+1:szA(2)+szB(2));
C = pdsi(szA(1)+1:szA(1)+szC(1), 1:szC(2));
D = pdsi(szA(1)+1:szA(1)+szC(1), szA(2)+1:szA(2)+szB(2));
% get the parameter space
pvs = psinfo(tpG,'par');
% find its size, skip zero sizes, add 1 for the nominal model
for szpvs=1:size(pvs,1)
    if pvs(szpvs,3)==0
        break;
    end
end
% loop through the vertices (~parameters ish). n is the local vertex index
for n=1:szpvs-1 
    % loop through the grid 
    gridres = 10;
    delta = (pvs(n,3)-pvs(n,2))/gridres;
    for ipvs=1:gridres % i.e grid split into 10 parts, ipvs is the grid index for this vertex
        % increase grid point value
        gpd = pvs(n,2)+(delta*ipvs); 
        % extract model at the point. polydec converts to polytopic form
        gp = polydec(pv,gpd);
        pdm = psinfo(tpG,'eval',gp);
        % it psinfo return an ltisys, this extracts the A matrix from it
        A(:,:,ipvs) = pdm(1:szA(1),1:szA(2));
    end
    % If multiple vertices were supported A(:,:,:) could be repackaged here
    % into the form A(col,row,grid,vertex)
    % A2 = zeros(szA(1),szA(2),szpvs-1,gridres); % move this!
    % A2(:,:,n,:) = A(:,:,:);
end

%%
% 3. Apply weights

% Design weights
s = tf('s');
W1 = tf([1.3],[0.5 5.9]);
W2 = 0.4;

% apply using series and ss
A1 = A; % A will need to be a different size anyway
clear('A')
for i=1:size(A1,3)
    Gs = ss(A1(:,:,i),B,C,D);
    Gs = series(series(Gs,W2),W1);
    % put back in the pv 'A' matrix
    [Ai,Bi,Ci,Di] = ssdata(Gs);
    % recreate
    if exist('A','var')==0
        A = zeros(size(Ai,1), size(Ai,2), size(A1,3));
    end
    A(:,:,i) = Ai;
end
% B,C,D also different, though not fundamentally, just aug'd 
B = Bi;
C = Ci;
D = Di;
% 

%%
% 4. Solve coprime factorisation for P, Z

% R
[rowd,cold]=size(D);
Rt=eye(rowd)+D*D';
R=eye(cold)+D'*D;

% initialise the LMIs
setlmis([]);
[rowp,colp]=size(B);
P=lmivar(1,[rowp 1]);
Z=lmivar(1,[rowp 1]);

% Loop through all varying plants in A
% one LMI for each Ai matrix ie an A model exists for each vertex
% B, C and D have to be fixed, so no looping takes place
for i=1:size(A,3)
    % first lmi term
    jj=newlmi;
    lmiterm([jj,1,1,P], 1,A(:,:,i),'s');
    lmiterm([jj,1,1,0], -C'*C);
    lmiterm([jj,2,1,P],B',1);
    lmiterm([jj,2,1,0],-D'*C);
    lmiterm([jj,2,2,0], -R);
end
% second lmi term
jj=newlmi;
lmiterm([-jj 1 1 Z],1,1);
lmiterm([-jj 2 1 0],1);
lmiterm([-jj 2 2 P],1,1);
lmisys=getlmis;
c = mat2dec(lmisys,zeros(rowp), eye(rowp));
options = [1e-2 200 1e5 10 0]
[copt,xopt] = mincx(lmisys,c,options)

%%
% 5. Form the coprime factorisation for each grid point using the solution
% P and Z

Popt=dec2mat(lmisys,xopt,1);
Zopt=dec2mat(lmisys,xopt,2);
% this is from Prempain LPV control P2
% Note that the coprime factorisation solution Popt enters here!
LL=-(B*D'+inv(Popt)*C')*inv(Rt); 
[rowc,colc]=size(C(:,:,1));

%---contruct the coprime factorisation design structure
% This is where the coprime plant is produced, so it is analagous to G in a
% ss config. The Prempain formulation is similar but different. 
for i=1:size(A,3)
    Acop=A(:,:,i);
    Bcop=[-LL*sqrt(Rt) B(:,:)];
    c1=[zeros(size(C));C(:,:)];
    c2=C(:,:);
    Ccop=[c1;c2];
    d11 = [zeros(size(D)); sqrt(Rt)];
    d12 = [eye(size(D)); sqrt(Rt)];
    d21 = sqrt(Rt);
    d22 = D(:,:);
    Dcop=[d11 d12; d21 d22];
    ss2 = pck(Acop,Bcop,Ccop,Dcop);
    Gcop(:,:,i)=pck(Acop,Bcop,Ccop,Dcop);		% CPF LPV Model
    Gcop2(:,:,i) = ss(Acop,Bcop,Ccop,Dcop);     % try to use modern notation
end

%% 
% 6. Synthesise a H-infinity controller using standard tools
% try to run hinf control on this

r = [3 3];

for i=1:size(A,3)
    P = Gcop(:,:,i);
    [gopt1,K1] = hinflmi(P,r); % testing only
    %[gopt,K,X1,X2,Y1,Y2] = hinflmi(P,r) % full
    
    
    [Katmp, Kbtmp, Kctmp, Kdtmp] = unpck(K1);
    % controller needs to be augmented with the weights in reverse order
    %Ks = series(series(K1,W2),W2);
    K2 = ss(Katmp,Kbtmp,Kctmp,Kdtmp);
    K2 = series(series(K2,W1),W2);
    [Katmp,Kbtmp,Kctmp,Kdtmp] = ssdata(K2);
    % array of ccontrollers
    K(:,:,1) = K2;
    gopt(i) = gopt1;
    % array of controller components, easier to extract!
    Ka(:,:,i) = Katmp;
    Kb(:,:,i) = Kbtmp;
    Kc(:,:,i) = Kctmp;
    Kd(:,:,i) = Kdtmp;
end

% compare gamma results
plot(1:size(A,3),gopt(:))
ylabel('Gamma Optimality')
xlabel('Param Index')

% nominal step response
G2 = ss(A(:,:,1),B,C,D);
K2 = ss(Ka(:,:,1), Kb(:,:,1), Kc(:,:,1), Kd(:,:,1))
CL = feedback(G2,-K2);
step(CL)

% loop functions
loops = loopsens(G2,-K2);
sigma(loops.So,loops.To,loops.CSo,{1e-2 1e3});
legend('S','T','KS');
title('Loop Shapes: S, T and KS')
legend('S','T','KS')

% controller gain
sigma(K2,{1e-2 1e3})
title('Controller Gain')

% pole zero plot - all frequencies
pzp = pzplot(K2);

% pole zero plot - lower frequencies
pzp = pzplot(K2);
pzpopt = getoptions(pzp);
pzpopt.XLim = {[-400 0]};
setoptions(pzp,pzpopt);
title('Controller Pole Zero Plot')

% multiple responses
for i=1:size(A,3)
    G2 = ss(A(:,:,i),B,C,D);
    K2 = ss(Ka(:,:,i), Kb(:,:,i), Kc(:,:,i), Kd(:,:,i))
    CL = feedback(G2,-K2);
    
    % for step response
    step(CL)
    
    % loop functions
%     loops = loopsens(G2,-K2);
%     sigma(loops.So,loops.To,loops.CSo);
%     legend('S','T','KS');
%     title('Loop Shapes: S, T and KS')
    
    % for controller gain
%     sigma(K2)
%     title('Controller Gain')
    
    
    hold on
end
hold off
legend('1','2','3','4','5','6','7','8','9','10')


%%
% 7. Mixed Synthesis
% There are likely to be HF poles placed by the algorithm so to counter
% this use pole placement constraints to limit HF poles. Can also place a
% bound on H2 aka LQR 
%[gopt,h2opt,K,R,S] = hinfmix(P,r,obj,region,dkbnd,tol)
r = [0 3 3]; % dim h2, y, u
obj = [0 0 1 0]; % h8 up bnd, h2 up bdn, weighting of h8, weighting of h2 
%region1 = lmireg;
% disc, centre at 0, radius 200
region1 = 1.0e+02 * ...
  [-2.0000 + 0.0200i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i  -2.0000 + 0.0000i   0.0100 + 0.0000i   0.0000 + 0.0000i];
[gopt,h2opt,K,R,S] = hinfmix(P,r,obj,region1)

% nominal step response
G2 = ss(A(:,:,1),B,C,D);
[Katmp, Kbtmp, Kctmp, Kdtmp] = unpck(K);
K2 = ss(Katmp,Kbtmp,Kctmp,Kdtmp);
K2 = series(series(K2,W1),W2);
CL = feedback(G2,-K2);
step(CL)

% loop functions
loops = loopsens(G2,-K2);
sigma(loops.So,loops.To,loops.CSo,{1e-2 1e3});
legend('S','T','KS');
title('Loop Shapes: S, T and KS')
legend('S','T','KS')

% controller gain
sigma(K2,{1e-2 1e3})
title('Controller Gain')

% pole zero plot - all frequencies
pzplot(K2);


% synth for all parameters
for i=1:size(A,3)
    P = Gcop(:,:,i);
    %region2 = lmireg;
    region1 = 1.0e+02 * ...
  [-2.0000 + 0.0200i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i  -2.0000 + 0.0000i   0.0100 + 0.0000i   0.0000 + 0.0000i];
    [gopt1,h2opt,K,R,S] = hinfmix(P,r,obj,region1);

    [Katmp, Kbtmp, Kctmp, Kdtmp] = unpck(K);
    % controller needs to be augmented with the weights in reverse order
    %Ks = series(series(K1,W2),W2);
    K2 = ss(Katmp,Kbtmp,Kctmp,Kdtmp);
    K2 = series(series(K2,W1),W2);
    [Katmp,Kbtmp,Kctmp,Kdtmp] = ssdata(K2);
    % array of ccontrollers
    gopt(i) = gopt1;
    % array of controller components, easier to extract!
    Ka(:,:,i) = Katmp;
    Kb(:,:,i) = Kbtmp;
    Kc(:,:,i) = Kctmp;
    Kd(:,:,i) = Kdtmp;
end

% compare gamma results
plot(1:size(A,3),gopt(:))
ylabel('Gamma Optimality')
xlabel('Param Index')

% multiple responses
for i=1:size(A,3)
    G2 = ss(A(:,:,i),B,C,D);
    K2 = ss(Ka(:,:,i), Kb(:,:,i), Kc(:,:,i), Kd(:,:,i));
    CL = feedback(G2,-K2);
    
    % for step response
    %step(CL)
    
    % loop functions
%     loops = loopsens(G2,-K2);
%     sigma(loops.So,loops.To,loops.CSo);
%     legend('S','T','KS');
%     title('Loop Shapes: S, T and KS')
    
    % for controller gain
    sigma(K2)
    title('Controller Gain')
    
    
    hold on
end
hold off
legend('1','2','3','4','5','6','7','8','9','10')
