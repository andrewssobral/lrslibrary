% Serhat Selcuk Bucak, bucakser@msu.edu

function ssbinitial(K,R)

     % blok boyu, þimdilik 16x16'lýk bloklar kullanýyoruz

mkdir('initials');
fname = ['initials/Initials' num2str(K) 'x' num2str(R) '.txt'];
fid=fopen(fname,'w');
for i = 1:K
    for j = 1:R
        %WeightVectors = randn(i,(2*P*R1 + K)*R2);
        fprintf(fid,'%8.4f\t', rand*255);
    end
    fprintf(fid,'\n');
end
fclose(fid);
