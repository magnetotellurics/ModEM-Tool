function [status] = Test3D(todo,arg1,arg2,arg3,arg4,arg5,arg6)

%   figure out if this is a serial or parallel job
nProc = getenv('NMPIPROC');
if ~isempty(nProc)
    if str2num(nProc) > 2
        SoP = 'MPI';
        MPIRun = ['mpirun -n ' nProc];
    else
        SoP = 'Serial';
    end
else
    SoP  = 'Serial';
end

if nargin == 3
    args = [arg1 ' ' arg2];
elseif nargin == 4
    args = [arg1 ' ' arg2 ' ' arg3];
elseif nargin == 5
    args = [arg1 ' ' arg2 ' ' arg3 ' ' arg4];
elseif nargin == 6
    args = [arg1 ' ' arg2 ' ' arg3 ' ' arg4 ' ' arg5];
elseif nargin == 7
    args = [arg1 ' ' arg2 ' ' arg3 ' ' arg4 ' ' arg5 ' ' arg6];
else
    error('Wrong number of arguments to Test3D');
end

switch SoP
    case 'Serial'
        switch todo
            case 'READ_WRITE'
                disp(['Mod3DMT -R ' args]);
                status = system(['Mod3DMT -R ' args]);
            case 'FORWARD'
                disp(['Mod3DMT -F ' args]);
                status = system(['Mod3DMT -F ' args]);
            case 'COMPUTE_J'
                disp(['Mod3DMT -J ' args]);
                status = system(['Mod3DMT -J ' args]);
            case 'MULT_BY_J'
                disp(['Mod3DMT -M ' args]);
                status = system(['Mod3DMT -M ' args]);
            case 'MULT_BY_J_T'
                disp(['Mod3DMT -T ' args]);
                status = system(['Mod3DMT -T ' args]);
            case 'MULT_BY_J_T_Mtx'
                disp(['Mod3DMT -x ' args]);
                status = system(['Mod3DMT -x ' args]);
            case 'INVERSE_NLCG'
                disp(['Mod3DMT -I NLCG ' args]);
                status = system(['Mod3DMT -I NLCG ' args]);
            otherwise
                error('This job is not implemented');
        end
    case 'MPI'     
        switch todo
            case 'READ_WRITE'
                disp([MPIRun ' Mod3DMT -R ' args]);
                status = system([MPIrun ' Mod3DMT_MPI -R ' args]);
            case 'FORWARD'
                disp([MPIRun ' Mod3DMT_MPI -F ' args]);
                status = system([MPIRun ' Mod3DMT_MPI -F ' args]);
            case 'COMPUTE_J'
                disp([MPIRun ' Mod3DMT_MPI -J ' args]);
                status = system([MPIRun ' Mod3DMT_MPI -J ' args]);
            case 'MULT_BY_J'
                disp([MPIRun ' Mod3DMT_MPI -M ' args]);
                status = system([MPIRun ' Mod3DMT_MPI -M ' args]);
            case 'MULT_BY_J_T'
                disp([MPIRun ' Mod3DMT_MPI -T ' args]);
                status = system([MPIRun ' Mod3DMT_MPI -T ' args]);
             case 'MULT_BY_J_T_Mtx'
                disp([MPIRun 'Mod3DMT_MPI -x ' args]);
                status = system([MPIRun ' Mod3DMT_MPI -x ' args]);
            case 'INVERSE_NLCG'
                disp([MPIRun 'm Mod3DMT_MPI -I NLCG ' args]);
                status = system([MPIRun ' Mod3DMT_MPI -I NLCG ' args]);
            otherwise
                error('This job is not implemented');
        end
end