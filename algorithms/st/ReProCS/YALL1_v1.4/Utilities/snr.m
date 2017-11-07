function SNR = snr(sig, ref)
% Compute Signal-to-Noise Ratio for signal 
%         relative to a reference
%
% Usage:
%       SNR = snr(sig, ref)  -- 1st time or
%       SNR = snr(sig)       -- afterwards
%
% Input:
%       sig         signal
%       ref         Reference
%  
% Output:
%       SNR           SNR value

persistent ref_sav ref_var

if nargin == 2
    ref_sav = ref;
    ref_var = var(ref(:),1);
end

mse = mean((ref_sav(:)-sig(:)).^2);
SNR = 10*log10(ref_var/mse);
