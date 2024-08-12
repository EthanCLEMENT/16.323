% Transfer function
s = tf('s');
G = 200/(0.025*s^3+1.0025*s^2+10.1*s+1);
numerator = [200];
denominator = [0.025 1.0025 10.1 1]; 

% State-space
[A,B,C,D] = tf2ss(numerator, denominator);
sys = ss(A,B,C,D);

% Bode plot of G
figure
bode(G);

M = 1.5;
wB = 10;
a=1e-4;

% Weighing functions
[Aw,Bw,Cw,Dw]=tf2ss([1/M wB],[1 wB*a]);
W1=ss(Aw,Bw,Cw,Dw);

W2 = 1;
W3 = [];

Paug = [A [0 0; 0 0; 0 0] B;-Bw*C Aw Bw 0;-Dw*C Cw Dw 0; 0 0 0 0 0 W2; -C 0 1 0];


% Augmented Matrix
P = augtf(G,W1,W2,W3);

test = [0 W1;0 W2 * G; W3*eye(size(G)) W3*G; eye(size(G)) G];
new_sys = ss(test);

Paug2= pck(test.A,test.B,test.C,test.D);

% H-infinity synthesis controller 
[K,CL,gamma] = hinfsyn(Paug2,1,1);


% Bisectioned algorithm to compute H-infinity norm
gama_l = 0.1;
gama_u = 100;
tol = 1e-6;


for i = 1:50
    if (gama_u-gama_l)/gama_l <= tol , break , end
    
    gama = (gama_u+gama_l)/2 ;
    
    H = [A (1/gama^2)*B*B'; -C'*C -A'];
    eigH = eig(H)' ;


    if min( abs(real(eigH)) ) < 1.0e-5 , gama_l = gama ;
    else gama_u = gama ;
end

end
i,gama,H
disp(hinfnorm(G));
Sensitivity_function = 1/(1+G*K);
robustness = hinfnorm(Sensitivity_function*W1);

% x = Ax  + B1w  + B2u 
% z = C1x + D11w + D12u
% y = C2x + D21w + D22u


% Assume Paug has the form:
% Paug = [A, B1, B2;
%         C1, D11, D12;
%         C2, D21, D22];



