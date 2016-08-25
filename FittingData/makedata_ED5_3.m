Cb_part1 = zeros(11,11,11);%Qb/Qc,AbAc,AsAc
Cs_part1 = zeros(11,11,11);%Qs/Qc,AbAc,AsAc
Cb_part2 = zeros(11,11,11);%Qb/Qc,AbAc,AsAc
Cs_part2 = zeros(11,11,11);%Qs/Qc,AbAc,AsAc

FileName='ED5_3.txt';
for ii=1:11
    Cb_part1(:,:,ii) = Get3DTable('ED5_3.txt', (ii-1)*11+6,(ii)*11+5)';
    Cs_part1(:,:,ii) = Get3DTable('ED5_3.txt', (ii-1)*11+133,(ii)*11+132)';
    Cb_part2(:,:,ii) = Get3DTable('ED5_3.txt', (ii-1)*11+260,(ii)*11+259)';
    Cs_part2(:,:,ii) = Get3DTable('ED5_3.txt', (ii-1)*11+387,(ii)*11+386)';
end

QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
QsQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AsAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];

[N1,N2,N3] = size(Cb_part1);
Cb_part1 = Cb_part1.*repmat(reshape(QbQc',N1,1,1),1,N2,N3)./repmat(reshape(AbAc,1,N2,1),N1,1,N3).^2;
[N1,N2,N3] = size(Cb_part2);
Cb_part2 = Cb_part2.*repmat(reshape(QbQc',N1,1,1),1,N2,N3)./repmat(reshape(AbAc,1,N2,1),N1,1,N3).^2;

[N1,N2,N3] = size(Cs_part1);
Cs_part1 = Cs_part1.*repmat(reshape(QsQc',N1,1,1),1,N2,N3)./repmat(reshape(AsAc,1,1,N3),N1,N2,1).^2;
[N1,N2,N3] = size(Cs_part2);
Cs_part2 = Cs_part2.*repmat(reshape(QsQc',N1,1,1),1,N2,N3)./repmat(reshape(AsAc,1,1,N3),N1,N2,1).^2;

save('ED5_3.mat','Cb_part1','Cb_part2','Cs_part1','Cs_part2','QbQc','QsQc','AsAc','AbAc');