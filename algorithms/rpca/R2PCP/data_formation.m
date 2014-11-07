switch input_sig
    case 1          % 1: synthetic data
        if load_sig
            load matlab
        else
        A_0=randn(n1,data.r)*randn(data.r,n2);  % true low-rank matrix A_0
%         A_max=max(abs(A_0(:)));
        A_max=sqrt(max(n1,n2));
        B_0=((rand(n1,n2)<.5)-.5)*2*A_max;
        B_0(randperm(n1*n2,n1*n2-data.s))=0;    % true sparse matrix B_0
%         X_0=A_0+B_0; Z=X_0;
        Z=A_0+B_0;
%         Z=randn(size(Z));
        end
        
    case 2          % 2: airport video
        addpath('./video_airport')
        nFrame=200;
        iFrame=1000;
        pathname=['airport',num2str(iFrame),'.bmp'];
        Z1=imread(pathname);
        Z1=double(uint8(round(sum(Z1,3)/3)))/255;
        Z=zeros(size(Z1,1),size(Z1,2),nFrame);
        Z(:,:,1)=Z1;
        for j=1:nFrame-1
            pathname=['airport',num2str(iFrame+j),'.bmp'];
            Z1=imread(pathname);
            Z1=double(uint8(round(sum(Z1,3)/3)))/255;
            Z(:,:,j+1)=Z1;
        end
        
        [z1,z2,z3]=size(Z); 
        n1=z1*z2; n2=z3;
        Z=reshape(Z,n1,n2);
        
        clear nFrame iFrame pathname Z1
        
    case 3          % 3: lobby video
        addpath('./video_lobby')
        nFrame=400;
        iFrame=1900;
        pathname=['SwitchLight',num2str(iFrame),'.bmp'];
        Z1=imread(pathname);
        Z1=double(uint8(round(sum(Z1,3)/3)))/255;
        Z=zeros(size(Z1,1),size(Z1,2),nFrame);
        Z(:,:,1)=Z1;
        for j=1:nFrame-1
            pathname=['SwitchLight',num2str(iFrame+j),'.bmp'];
            Z1=imread(pathname);
            Z1=double(uint8(round(sum(Z1,3)/3)))/255;
            Z(:,:,j+1)=Z1;
        end
        
        [z1,z2,z3]=size(Z); 
        n1=z1*z2; n2=z3;
        Z=reshape(Z,n1,n2);
        
        clear nFrame iFrame pathname Z1
end
        
% add noise
% noise=.0;
% noise=.05*std(data);
if ~load_sig && noise~=0
    Z=Z+noise*randn(size(Z));                       % Gaussian noise
end
