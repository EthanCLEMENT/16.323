set(0,'DefaultAxesFontName','arial')
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultTextFontName','arial')
set(0,'DefaultTextFontSize',20)
clear all
if ~exist('yprev')
    yprev=[1 1]';
    tprev=[0 1]';
    Sensprev=[1 1];
    fprev=[.1 100];
end
%Wu=1/1e9;
Wu=1;

% define plant
[Ag,Bg,Cg,Dg]=tf2ss(200,conv(conv([0.05 1],[0.05 1]),[10 1]));
Gol=ss(Ag,Bg,Cg,Dg);

% define sensitivity weight
M=1.5;wB=10;A=1e-4;
[Asw,Bsw,Csw,Dsw]=tf2ss([1/M wB],[1 wB*A]);
Ws=ss(Asw,Bsw,Csw,Dsw);

% form augmented P dynamics
n1=size(Ag,1);
n2=size(Asw,1);
A=[Ag zeros(n1,n2);-Bsw*Cg Asw];
Bw=[zeros(n1,1);Bsw];
Bu=[Bg;zeros(n2,1)];
Cz=[-Dsw*Cg Csw;zeros(1,n1+n2)];
Cy=[-Cg zeros(1,n2)];
Dzw=[Dsw;0];
Dzu=[0;Wu];
Dyw=[1];
Dyu=0;
P=pck(A,[Bw Bu],[Cz;Cy],[Dzw Dzu;Dyw Dyu]);

% Given matrices
B = [Bw, Bu];
m1 = size(Bw,2);
m2 = size(Bu,2);
gamma = 1.365;

Q = Cz'*Cz;
R = [-gamma^2*eye(m1) zeros(m1,m2) ; zeros(m2,m1) eye(m2)];
X = icare(A,[Bw,Bu],Q,R);
disp(X);

% Solve the ARE
[X, L1, G1] = care(A, B, Q, R);

C = [Cz' Cy'];
m1 = size(Cz',2);
m2 = size(Cy',2);
Q = Bw'*Bw;
R =[-gamma^2*eye(m1) zeros(m1,m2) ; zeros(m2,m1) eye(m2)];

% Solve the ARE
[Y, L2, G2] = care(A, C, Q, R);

disp(max(abs(eig(X*Y))) < gamma^2);
Z = inv((eye(4) - gamma^-2*X*Y));

HA = [A+(gamma^-2*Bw*Bw' - Bu*Bu')*X - Z*Y*Cy'*Cy];
HB = [Z*Y*Cy'];
HC = [-Bu'*X];
HD = 0;
disp([A+(gamma^-2*Bw*Bw' - Bu*Bu')*X - Z*Y*Cy'*Cy Z*Y*Cy';-Bu'*X  0]);
gama_l = 0.1;
gama_u = 100;
tol = 1e-6;
for i = 1:50
    if (gama_u-gama_l)/gama_l <= tol , break , end
    
    gama = (gama_u+gama_l)/2 ;
    
    H = [HA (1/gama^2)*HB*HB'; -HC'*HC -HA'];
    eigH = eig(H)' ;


    if min( abs(real(eigH)) ) < 1.0e-5 , gama_l = gama ;
    else gama_u = gama ;
end

end
i,gama,H

% call hinf to find Gc (mu toolbox)
diary hinf1_diary
[Gc,G,gamma]=hinfsyn(P,1,1,0.1,20,.001);
diary off
[ac,bc,cc,dc]=unpck(Gc);
ev=max(real(eig(ac)/2/pi))
PP=ss(A,[Bw Bu],[Cz;Cy],[Dzw Dzu;Dyw Dyu]);
GGc=ss(ac,bc,cc,dc);
CLsys = feedback(PP,GGc,[2],[3],1);
[acl,bcl,ccl,dcl]=ssdata(CLsys);

% reduce closed-loop system so that it only has 1 input and 2 outputs
bcl=bcl(:,1);ccl=ccl([1 2],:);dcl=dcl([1 2],1);
CLsys=ss(acl,bcl,ccl,dcl);
f=logspace(-1,2,400);
Pcl=freqresp(CLsys,f);
CLWS=squeeze(Pcl(1,1,:)); % closed loop weighted sens
WS=freqresp(Ws,f); % sens weight
SensW=squeeze(WS(1,1,:));
Sens=CLWS./SensW; % divide out weight to get closed-loop sens

figure(1);clf
loglog(f,abs(Sens),'b-','LineWidth',2)
hold on
loglog(f,abs(1./SensW),'m--','LineWidth',2)
loglog(f,abs(CLWS),'r-.','LineWidth',2)
loglog(fprev,abs(Sensprev),'r.')
legend('S','1/W_s','W_sS','Location','SouthEast')
hold off
xlabel('Freq (rad/sec)')
ylabel('Magitude')
grid


na=size(Ag,1);
nac=size(ac,1);
Acl=[Ag Bg*cc;-bc*Cg ac];Bcl=[zeros(na,1);bc];Ccl=[Cg zeros(1,nac)];Dcl=0;
Gcl=ss(Acl,Bcl,Ccl,Dcl);
[y,t]=step(Gcl,1);

figure(2);clf
plot(t,y,'LineWidth',2)
hold on;plot(tprev,yprev,'r--','LineWidth',2);hold off
xlabel('Time sec')
ylabel('Step response')

yprev=y;
tprev=t;
Sensprev=Sens;
fprev=f;


