function stop = outfun(x, optimValues, state) % output function

stop = false;
grad = optimValues.gradient;

fileName = 'grad_output.txt';
[fid, message] = fopen(fileName,'a');

switch state
    
    case 'init'
        
        if fid==-1
            error('globaloptim:psoutputfile:fileError','Error trying to write to %s:\n%s',fileName,message)
        end
        
    case 'iter'
        
        Iter = optimValues.iteration;
        % Write to the file
        fprintf(fid,'Iter: %d \t',Iter);
        for i=1:length(grad)
            fprintf(fid,'%5.4f \t', grad(i));
        end
        fprintf(fid,'\n');
        
    case 'done'
        
        Iter = optimValues.iteration;
        % Write to the file
        fprintf(fid,'Iter: %d \t',Iter);
        for i=1:length(grad)
            fprintf(fid,'%5.4f \t', grad(i));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\nOptimization terminated.\n');
        fclose(fid);
        
end

