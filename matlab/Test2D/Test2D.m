function [status] = Test2D(todo,arg1,arg2,arg3,arg4,arg5,arg6)

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
    error('Wrong number of arguments to Test2D');
end

switch todo
    case 'READ_WRITE'
        disp(['Mod2DMT -R ' args]);
        status = system(['Mod2DMT -R ' args]);
    case 'FORWARD'
        disp(['Mod2DMT -F ' args]);
        status = system(['Mod2DMT -F ' args]);
    case 'COMPUTE_J'
        disp(['Mod2DMT -J ' args]);
        status = system(['Mod2DMT -J ' args]);
    case 'MULT_BY_J'
        disp(['Mod2DMT -M ' args]);
        status = system(['Mod2DMT -M ' args]);
    case 'MULT_BY_J_T'
        disp(['Mod2DMT -T ' args]);
        status = system(['Mod2DMT -T ' args]);
    case 'INVERSE_NLCG'
        disp(['Mod2DMT -I NLCG ' args]);
        status = system(['Mod2DMT -I NLCG ' args]);
    otherwise
        error('This job is not implemented');
end
