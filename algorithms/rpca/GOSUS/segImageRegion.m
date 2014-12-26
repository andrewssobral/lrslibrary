function I_s = segImageRegion(I,S)

I_s = I;
I_s1 = I(:,:,1);
I_s2 = I(:,:,2); 
I_s3 = I(:,:,3); 

for i=1:max(S(:))
   
    mask = S==i;
    
%     I_s(mask) = mean(I(mask));
    I_s1(mask) = mean(I_s1(mask));
    I_s2(mask) = mean(I_s2(mask));
    I_s3(mask) = mean(I_s3(mask));

    I_s(:,:,1) = I_s1;
    I_s(:,:,2) = I_s2;
    I_s(:,:,3) = I_s3;
end

end
