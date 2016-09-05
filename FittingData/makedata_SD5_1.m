Cb = zeros(11,11);%Qb/Qc,AbAc
Cs = zeros(11,11);%Qs/Qc,AsAc

FileName='SD5_1.txt';

Cb = Get2DTable(FileName, 6,16)';
Cs = Get2DTable(FileName, 23,33)';

QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
QsQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AsAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];

[N1,N2] = size(Cb);
Cb = Cb.*repmat(reshape(QbQc',N1,1),1,N2)./repmat(reshape(AbAc,1,N2),N1,1).^2;
[N1,N2] = size(Cs);
Cs = Cs.*repmat(reshape(QsQc',N1,1),1,N2)./repmat(reshape(AsAc,1,N2),N1,1).^2;

save('SD5_1.mat','Cb','Cs','QbQc','QsQc','AsAc','AbAc');
