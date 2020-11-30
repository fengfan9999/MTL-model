function write_comsol_grid(filename,data,x,y,z)

% if (nargin<4)
%     error('Not enough input arguments. Syntax is \n\n write_comsol_grid(filename,data,x,y)\n\n or \n\n write_comsol_grid(filename,data,x,y,z)',1)
% end

if nargin==2
    f=fopen(filename,'a');
else
    f=fopen(filename,'w');
end

if nargin>2
    
    disp('Writing Grid...');
    fprintf(f,'%%Grid \r\n');
    fprintf(f, [num2str(x(:)') '\r\n']);
    fprintf(f, [num2str(y(:)') '\r\n']);
    if (nargin>4)
        fprintf(f, [num2str(z(:)') '\r\n']);
    end

    disp('Writing Data...');
else
    disp('Appending Data');
end


data=permute(data,[2 1 3]);


fprintf(f,'%%Data \r\n');
fprintf(f, [num2str(data(:)') '\r\n']);

fclose(f);

