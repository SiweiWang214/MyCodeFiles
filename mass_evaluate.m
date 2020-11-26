function f=mass_evaluate(f, M, V)
%evaluate objectives for multiple chromosome in parallel by submitting jobs
%to LSF. Modify function submit_job, check_job and get_job_data to use.
%
%f,  N x (V+M), chromosome
%M,  integer, number of objectives
%V,  integer, number of variables
% created 10/12/2011, Xiaobiao Huang
%

N = size(f,1);
%if M~=1
%    error('there can be only 1 objective');
%end
global g_cnt g_data

parfor i=1:N
    x = f(i,1:V);
    
    obj1 = func_obj(x);
    if mod(g_cnt,10)==0
        [g_cnt, x(:)', obj1];
    end
    
    res_tmp(i,:) = obj1;
%    f(i,V + 1: M + V) = obj1;
end

f(:,V + 1: M + V) = res_tmp;


