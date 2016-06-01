function  [rows, cols, entries, rep, field, symm] = mminfo(filename)
%
%  function  [rows, cols, entries, rep, field, symmetry] = mminfo(filename)
%
%      Reads the contents of the Matrix Market file 'filename'
%      and extracts size and storage information.
%
%      In the case of coordinate matrices, entries refers to the
%      number of coordinate entries stored in the file.  The number
%      of non-zero entries in the final matrix cannot be determined
%      until the data is read (and symmetrized, if necessary).
%
%      In the case of array matrices, entries is the product
%      rows*cols, regardless of whether symmetry was used to
%      store the matrix efficiently.
%
%

mmfile = fopen(filename,'r');
if ( mmfile == -1 )
 disp(filename);
 error('File not found');
end;

header = fgets(mmfile);
if (header == -1 )
  error('Empty file.')
end
      
% NOTE: If using a version of Matlab for which strtok is not
%       defined, substitute 'gettok' for 'strtok' in the 
%       following lines, and download gettok.m from the
%       Matrix Market site.    
[head0,header]   = strtok(header);  % see note above
[head1,header]   = strtok(header);
[rep,header]     = strtok(header);
[field,header]   = strtok(header);
[symm,header]    = strtok(header);
head1 = lower(head1);
rep   = lower(rep);
field = lower(field);
symm  = lower(symm);
if ( length(symm) == 0 )
   disp('Not enough words in header line.') 
   disp('Recognized format: ')
   disp('%%MatrixMarket matrix representation field symmetry')
   error('Check header line.')
end
if ( ~ strcmp(head0,'%%MatrixMarket') )
   error('Not a valid MatrixMarket header.')
end
if (  ~ strcmp(head1,'matrix') )
   disp(['This seems to be a MatrixMarket ',head1,' file.']);
   disp('This function only knows how to read MatrixMarket matrix files.');
   disp('  ');
   error('  ');
end

% Read through comments, ignoring them

commentline = fgets(mmfile);
while length(commentline) > 0 & commentline(1) == '%',
  commentline = fgets(mmfile);
end

% Read size information, then branch according to
% sparse or dense format

if ( strcmp(rep,'coordinate')) %  read matrix given in sparse 
                              %  coordinate matrix format

  [sizeinfo,count] = sscanf(commentline,'%d%d%d');
  while ( count == 0 )
     commentline =  fgets(mmfile);
     if (commentline == -1 )
       error('End-of-file reached before size information was found.')
     end
     [sizeinfo,count] = sscanf(commentline,'%d%d%d');
     if ( count > 0 & count ~= 3 )
       error('Invalid size specification line.')
     end
  end
  rows = sizeinfo(1);
  cols = sizeinfo(2);
  entries = sizeinfo(3);

elseif ( strcmp(rep,'array') ) %  read matrix given in dense 
                               %  array (column major) format

  [sizeinfo,count] = sscanf(commentline,'%d%d');
  while ( count == 0 )
     commentline =  fgets(mmfile);
     if (commentline == -1 )
       error('End-of-file reached before size information was found.')
     end
     [sizeinfo,count] = sscanf(commentline,'%d%d');
     if ( count > 0 & count ~= 2 )
       error('Invalid size specification line.')
     end
  end
  rows = sizeinfo(1);
  cols = sizeinfo(2);
  entries = rows*cols;
end

fclose(mmfile);
% Done.

