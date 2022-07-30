function u_out = RemoveWhiteSpace(u_in, varargin)
%  February 2nd, 2012, By Reza Farrahi Moghaddam, Synchromedia Lab, ETS, Montreal, Canada
%
% RemoveWhiteSpace function removes white spaces around an image.
%
% Syntax:
%	1. For an image: u_out = RemoveWhiteSpace(u_in)
%
%	2. For an image file, to write the result on the same file: RemoveWhiteSpace([], 'file', input_filename)
%
%	3. For an image file, to make a new output file: RemoveWhiteSpace([], 'file', input_filename, 'output', output_filename)
%

% get the arguments
[it_is_a_file_flag, input_filename, output_filename] = check_the_argin_infile(nargin, varargin{:});

%
if (it_is_a_file_flag)
	[u_in, map] = imread(input_filename);
	if (numel(map) ~= 0)
		u_in = ind2rgb(u_in, map);
	end
end
u_in = mat2gray(u_in);
[xm ym zm] = size(u_in);

%
if (zm == 3)
	u_gray = rgb2gray(u_in);
else
	u_gray = mean(u_in, 3);
end

%
u_white_mask = u_gray > 0.99;
u_white_mask_hori = reshape(sum(u_white_mask, 1), [], 1);
u_white_mask_vert = sum(u_white_mask, 2);
%
u_white_mask_hori(u_white_mask_hori < xm / 2) = 0;
u_white_mask_vert(u_white_mask_vert < ym / 2) = 0;
u_white_mask_hori_diff = diff(u_white_mask_hori);
u_white_mask_vert_diff = diff(u_white_mask_vert);
[~, boundingbox_hori] = findpeaks(abs(u_white_mask_hori_diff));
[~, boundingbox_vert] = findpeaks(abs(diff(u_white_mask_vert)));
if (numel(boundingbox_hori) == 0)
	boundingbox_hori = [0 ym]; 
elseif (numel(boundingbox_hori) == 1)
	if (boundingbox_hori > ym / 2)
		boundingbox_hori = [0, boundingbox_hori];		
		
	else
		boundingbox_hori = [boundingbox_hori, ym];
	end
else
	boundingbox_hori = boundingbox_hori([1, end]);
	boundingbox_hori = boundingbox_hori .* [- u_white_mask_hori_diff(boundingbox_hori(1)) > 0; u_white_mask_hori_diff(boundingbox_hori(2)) > 0];
end
if (numel(boundingbox_vert) == 0)
	boundingbox_vert = [0 xm]; 
elseif (numel(boundingbox_vert) == 1)
	if (boundingbox_vert > xm / 2)
		boundingbox_vert = [0, boundingbox_vert];		
	else
		boundingbox_vert = [boundingbox_vert, xm];
	end	
else
	boundingbox_vert = boundingbox_vert([1, end]);
	boundingbox_vert = boundingbox_vert .* [- u_white_mask_vert_diff(boundingbox_vert(1)) > 0; u_white_mask_vert_diff(boundingbox_vert(2)) > 0];
end
boundingbox_hori(1) = boundingbox_hori(1) + 1;
boundingbox_vert(1) = boundingbox_vert(1) + 1;
if (boundingbox_hori(2) == 0)
	boundingbox_hori(2) = ym;
end
if (boundingbox_vert(2) == 0)
	boundingbox_vert(2) = xm;
end
%
u_out = u_in(boundingbox_vert(1) : boundingbox_vert(2), boundingbox_hori(1) : boundingbox_hori(2), :);

%
if (it_is_a_file_flag)
	imwrite(u_out, output_filename);
end

end


function [it_is_a_file_flag, input_filename, output_filename] = check_the_argin_infile(nargin, varargin)
% 120202: Reza

% %
it_is_a_file_flag = false;
input_filename = '';
output_filename = '';

%
default_fields = {'file', 'output'};
for temp_label = 1 : 2 : (nargin - 1)
    parameter_name = varargin{temp_label};
    parameter_val = varargin{temp_label + 1};
	matched_field = find(strcmpi(parameter_name, default_fields)); 
    if isempty(matched_field)
        error('Error: Unknown argument: %s.\n', parameter_name);
    else
        switch(matched_field)
            case 1  % input_filename
                input_filename = parameter_val;
                output_filename = parameter_val;				
				it_is_a_file_flag = true;
            case 2  % output_filename
                output_filename = parameter_val;
				it_is_a_file_flag = true;				
        end
    end
end

%
if (numel(output_filename) > 0)&&(numel(input_filename) == 0)
	input_filename = output_filename;
end

end


%{
u_in = mat2gray(imread('test.png'));
subfigure(2, 2, [2 2]), imshow(u_in, 'InitialMagnification','fit');
subfigure(2, 2, [1 1]), imshow(RemoveWhiteSpace(u_in), 'InitialMagnification','fit');
%}

%{
RemoveWhiteSpace([], 'file', 'test.png', 'output', 'test_out.png');
%}

%{
RemoveWhiteSpace([], 'file', 'test.png');
%}