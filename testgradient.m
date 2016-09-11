n = 10000;
dPdS_Jacobian = zeros(n,1);
Pdrop =  zeros(n,1);
dPdQ =  zeros(n,1);
dPdS_Cir =  zeros(n,6);
error =  zeros(n,1);
Q = rand(n,1);
L = 10*rand(n,1);
H = 0.5+0.5*rand(n,1);
W = 0.5+0.5*rand(n,1);
D = 0.5+0.5*rand(n,1);
rho = 1.204*(0.7+0.6*rand(n,1));
e = 0.09*(0.7+0.6*rand(n,1));
nu = 15.11e-6*(0.7+0.6*rand(n,1));
m = 2;
for i= 1:n
    dPdS_Jacobian(i) = DuctNetwork.Jacobian(@(x)DuctNetwork.RectangularDarcyWeisbach({'Pdrop','dPdQ','dPdS'},Q(i),[L(i),H(i),W(i),rho(i),e(i),nu(i)]+sparse(1,m,x,1,6)),0,1e-6*[L(i),H(i),W(i),rho(i),e(i),nu(i)]*sparse(m,1,1,6,1));
    [Pdrop(i),dPdQ(i),dPdS_Cir(i,:)] = DuctNetwork.RectangularDarcyWeisbach({'Pdrop','dPdQ','dPdS'},Q(i),[L(i),H(i),W(i),rho(i),e(i),nu(i)]);
    error(i)= ((dPdS_Cir(i,m)-dPdS_Jacobian(i))/dPdS_Cir(i,m));
end
[max_error,max_error_ID] = max(error);
i = max_error_ID;
[Pdrop_t,dPdQ_t,dPdS_Cir_t] = DuctNetwork.RectangularDarcyWeisbach({'Pdrop','dPdQ','dPdS'},Q(i),[L(i),H(i),W(i),rho(i),e(i),nu(i)]);
plot(error);