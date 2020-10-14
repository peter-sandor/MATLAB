function t_stat = student_t(a,v)
% calculate the Student's ratio for scalar 'a' and array 'v'
% (e.g. for testing whether the mean of the values in 'v' is equal to 'a')
N=length(v);
t_stat=sqrt(N/(N+1))*(a-mean(v))/std(v);
end