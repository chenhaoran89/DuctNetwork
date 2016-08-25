Cb = zeros(11,11);%Qb/Qc,AbAc
FileName='SD5_18.txt';
Cb = Get2DTable(FileName, 6,16)';
QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
[N1,N2] = size(Cb);
Cb = Cb.*repmat(reshape(QbQc',N1,1),1,N2)./repmat(reshape(AbAc,1,N2),N1,1).^2;
save('SD5_18.mat','Cb','QbQc','AbAc');

