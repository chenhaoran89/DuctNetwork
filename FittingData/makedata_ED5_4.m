Cb1 = zeros(11,11,11);%Qb/Qc,Ab1Ac,Ab2Ac
Cb2 = zeros(11,11,11);%Qs/Qc,Ab1Ac,Ab2Ac

FileName='ED5_4.txt';
for ii=1:11
    Cb1(:,:,ii) = Get3DTable('ED5_4.txt', (ii-1)*11+6,(ii)*11+5)';
    Cb2(:,:,ii) = Get3DTable('ED5_4.txt', (ii-1)*11+133,(ii)*11+132)';
end

QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
[N1,N2,N3] = size(Cb1);
Cb1 = Cb1.*repmat(reshape(QbQc',N1,1,1),1,N2,N3)./repmat(reshape(AbAc,1,N2,1),N1,1,N3).^2;
[N1,N2,N3] = size(Cb2);
Cb2 = Cb2.*repmat(reshape(QbQc',N1,1,1),1,N2,N3)./repmat(reshape(AbAc,1,1,N3),N1,N2,1).^2;

save('ED5_4.mat','Cb1','Cb2','QbQc','AbAc');