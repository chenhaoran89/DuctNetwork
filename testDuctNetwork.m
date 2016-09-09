duct = DuctNetwork();
duct.AddParameter({'Density(kg/m^3)','Roughness(mm)','Kinematic Viscosity(m^2/s)'},[1.204,0.15,15.11e-6]);
duct.AddBranch(0,1);
duct.AddFitting(1,'FanQuadratic',[150*(1+rand),0.3]);
duct.AddFitting(1,'DuctQuadratic',1000*[2+rand]);
duct.AddFitting(1,'CircularDarcyWeisbach',[2.2+rand,0.2]);
duct.AddFitting(1,'CircularDarcyWeisbach',[2.1+rand,0.2]);
duct.AddBranch(1,2);
duct.AddFitting(2,'DuctQuadratic',1000*[2+rand]);
duct.AddFitting(2,'CircularDarcyWeisbach',[2.2+rand,0.2]);
duct.AddFitting(2,'CircularDarcyWeisbach',[2.1+rand,0.2]);
duct.AddBranch(1,3);
duct.AddFitting(3,'DuctQuadratic',1000*[2+rand]);
duct.AddFitting(3,'CircularDarcyWeisbach',[2.2+rand,0.2]);
duct.AddFitting(3,'CircularDarcyWeisbach',[2.1+rand,0.2]);
duct.AddBranch(2,0);
duct.AddFitting(4,'DuctQuadratic',1000*[2+rand]);
duct.AddFitting(4,'CircularDarcyWeisbach',[2.2+rand,0.2]);
duct.AddFitting(4,'CircularDarcyWeisbach',[2.1+rand,0.2]);
duct.AddBranch(3,0);
duct.AddFitting(5,'DuctQuadratic',1000*[2+rand]);
duct.AddFitting(5,'CircularDarcyWeisbach',[2.1+rand,0.2]);
duct.AddBranch(2,duct.AddNode('Node 1'));
duct.AddFitting(6,'DuctQuadratic',1000*[2+rand]);
duct.AddFitting(6,'CircularDarcyWeisbach',[2.2+rand,0.2]);
duct.AddFitting(6,'CircularDarcyWeisbach',[2.1+rand,0.2]);
duct.AddBranch(3,duct.AddNode('Node 2'));
duct.AddFitting(7,'DuctQuadratic',1000*[2+rand]);
duct.AddFitting(7,'CircularDarcyWeisbach',[2.2+rand,0.2]);
duct.AddFitting(7,'CircularDarcyWeisbach',[2.1+rand,0.2]);
duct.AddBranch('Node 1','ATM');
duct.AddFitting(8,'DuctQuadratic',1000*[2+rand]);
duct.AddFitting(8,'CircularDarcyWeisbach',[2.2+rand,0.2]);
duct.AddFitting(8,'CircularDarcyWeisbach',[2.1+rand,0.2]);
duct.AddBranch('Node 2','GND');
duct.AddFitting(9,'DuctQuadratic',1000*[2+rand]);
duct.AddFitting(9,'CircularDarcyWeisbach',[2.2+rand,0.2]);
duct.AddFitting(9,'CircularDarcyWeisbach',[2.1+rand,0.2]);
duct.AddFitting([1,2,3],'CircularTJunction',[0.3,0.3,0.2]);

%%
N=100;
for ii=1:N
    disp(ii)
    duct.Sim();
end
SuccessRatio=2-duct.n_trail/N