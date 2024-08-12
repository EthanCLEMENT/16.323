A = [-1 3; -1 -6];
B = [0;1];
C = [2 -1];
D = [0];

sys = ss(A, B, C, D);

G = ss2tf(A,B,C,D);

integrandPc = @(t) expm(A'*t) * C' * C * expm(A*t);
ObservabilityGramian = integral(integrandPc, 0, inf, 'ArrayValued', true);

integrandPo = @(t) expm(A*t) * B * B' * expm(A'*t);
ControllabilityGramian = integral(integrandPo, 0, inf, 'ArrayValued', true);

Pc = gram(sys,'c');
Po = gram(sys,'o');

n2 = norm(sys,2);

H2 = trace(B'*Po*B);

ninf = hinfnorm(sys);

gama_l = 0.1;
gama_u = 5;
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



Ps = [1/(s+1) -1/(s^2) ; 1 -1/(s^2)];
new_sys = ss(Ps);
