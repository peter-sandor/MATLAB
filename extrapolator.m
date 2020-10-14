function varargout=extrapolator
temp=ginput(2);
temp=sortrows(temp,1);
S.userdata=temp;
S.func='y=A*x+B';
S.A=(temp(2,2)-temp(1,2))/(temp(2,1)-temp(1,1));
S.B=temp(2,2)-S.A*temp(2,1);
varargout{1}=S;
end