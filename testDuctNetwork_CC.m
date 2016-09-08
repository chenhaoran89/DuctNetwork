duct = DuctNetwork();
duct.AddParameter({'Density(kg/m^3)','Roughness(mm)','Dynamic Viscosity(m^2/s)'},[1.204,0.15,15.11e-6],[2,3,4]);
duct.AddBranch(0,1);
duct.AddFitting(1,'PressureSource',[20]);
duct.AddBranch(1,0);
duct.AddBranch(1,0);
duct.AddFitting([1,2,3],'SR5_13_Junction',[0.2,0.3,0.2,0.1,0.2,0.2]);
[X,Q,P]=duct.Sim();
%%
h=1e-6;
X0 = [0.9;0.22];
Branch = 3;
[dP,dPdQ,dPdS]=duct.BranchPressureDrop(Branch,X0,duct.S);
dP1 = duct.BranchPressureDrop(Branch,(X0+[h;0]),duct.S);
dP2 = duct.BranchPressureDrop(Branch,(X0+[0;h]),duct.S);
dP0dQ = [(dP1-dP)/h,(dP2-dP)/h];
dPdQ
dP0dQ
