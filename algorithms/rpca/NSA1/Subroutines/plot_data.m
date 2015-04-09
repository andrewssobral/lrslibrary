function plot_data(frame_array,D,X,S,n1,n2)

    T=size(frame_array,2);
    for k=1:T
        subplot('Position',[(k-1)/T 0.68 1/T 0.3]), imshow(reshape(D(:,frame_array(k)),n1,n2),[])
        subplot('Position',[(k-1)/T 0.35 1/T 0.3]), imshow(reshape(X(:,frame_array(k)),n1,n2),[])
        subplot('Position',[(k-1)/T 0.02 1/T 0.3]), imshow(reshape(S(:,frame_array(k)),n1,n2),[])
    end

% N=size(frame_array,2);
% for k=1:N
%     subplot(3,N,k), imshow(reshape(D(:,frame_array(k)),144,176),[])
%     subplot(3,N,k+N), imshow(reshape(X(:,frame_array(k)),144,176),[])
%     subplot(3,N,k+2*N), imshow(reshape(S(:,frame_array(k)),144,176),[])
% end