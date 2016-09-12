n = 1000;
dPdS_Jacobian = zeros(n,1);
Pdrop =  zeros(n,1);
dPdQ =  zeros(n,1);
dPdS_Rect =  zeros(n,6);
% dPdS_Cir =  zeros(n,5);
error =  zeros(n,1);
Q = rand(n,1);
L = 10*rand(n,1);
H = 0.5+0.5*rand(n,1);
W = 0.5+0.5*rand(n,1);
% D = 0.5+0.5*rand(n,1);
rho = 1.204*(0.7+0.6*rand(n,1));
e = 0.00009*(0.7+0.6*rand(n,1));
nu = 15.11e-6*(0.7+0.6*rand(n,1));
for i= 1:n
    dPdS_Jacobian(i) = DuctNetwork.Jacobian(@(x)DuctNetwork.RectangularDarcyWeisbachHaaland({'Pdrop','dPdQ','dPdS'},Q(i),[L(i),H(i),W(i),rho(i),e(i),x]),nu(i),1e-8);
    [Pdrop(i),dPdQ(i),dPdS_Rect(i,:)] = DuctNetwork.RectangularDarcyWeisbachHaaland({'Pdrop','dPdQ','dPdS'},Q(i),[L(i),H(i),W(i),rho(i),e(i),nu(i)]);
    error(i)= dPdS_Rect(i,6)/dPdS_Jacobian(i)-1;
    if abs(error(i))> 5e-4
        Dh = 1.3*H(i)^0.625*W(i)^0.625/(H(i)+W(i))^0.25;
        Re = abs(Q(i))*Dh/H(i)/W(i)/nu(i)
        lambda = 1/(1+exp(-(Re-3750)/250))
        A = (e(i)/3.7/Dh)^3.33
        B = (6.9/Re)^3
        [dPdS_Rect(i,6),dPdS_Jacobian(i),error(i)]
        
        %         [(e(i)/3.7/D(i))^3.33,(6.9/Re)^3]
        %         [Q(i),Q(i)/D(i)^2,4*abs(Q(i))/nu/pi/D(i),dPdS_Cir(i,4),dPdS_Jacobian(i)]
    end
end
% for i= 1:n
%     dPdS_Jacobian(i) = DuctNetwork.Jacobian(@(x)DuctNetwork.CircularDarcyWeisbachHaaland({'Pdrop','dPdQ','dPdS'},Q(i),[L(i),D(i),rho(i),x,nu(i)]),e(i),1e-6);
%     [Pdrop(i),dPdQ(i),dPdS_Cir(i,:)] = DuctNetwork.CircularDarcyWeisbachHaaland({'Pdrop','dPdQ','dPdS'},Q(i),[L(i),D(i),rho(i),e(i),nu(i)]);
%     error(i)= dPdS_Cir(i,4)/dPdS_Jacobian(i)-1;
%     if abs(error(i))> 1e-4
% %         Re= 4*abs(Q(i))/nu(i)/pi/D(i);
%         [dPdS_Cir(i,4),dPdS_Jacobian(i)]
% %         [(e(i)/3.7/D(i))^3.33,(6.9/Re)^3]
% %         [Q(i),Q(i)/D(i)^2,4*abs(Q(i))/nu/pi/D(i),dPdS_Cir(i,4),dPdS_Jacobian(i)]
%     end
% end
% [Pdrop_t,dPdQ_t,dPdS_Rect_t] = DuctNetwork.RectangularDarcyWeisbach({'Pdrop','dPdQ','dPdS'},Q(i),[L(i),H(i),W(i),rho(i),e(i),nu(i)]);
% [Pdrop_t,dPdQ_t,dPdS_Cir_t] = DuctNetwork.CircularDarcyWeisbachHaaland({'Pdrop','dPdQ','dPdS'},Q(i),[L(i),D(i),rho(i),e(i),nu(i)]);
plot(error);