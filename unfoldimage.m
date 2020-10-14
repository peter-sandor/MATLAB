function out = unfoldimage(in,quadrant)

% quadrant convention: top left: 1, top right: 2, bottom left: 3, bottom
% right: 4

switch quadrant
    case 1
        in = flipdim(flipdim(in,1),2);
    case 2
        in = flipdim(in,1);
    case 3
        in = flipdim(in,2);
    case 4
        % do nothing;
end
N=size(in);

M=2*N;
out=zeros(M);
out(N(1)+1:M(1),N(2)+1:M(2))=in;
out(1:N(1),N(2)+1:M(2))=flipdim(in,1);
out(N(1)+1:M(1),1:N(2))=flipdim(in,2);
out(1:N(1),1:N(2))=flipdim(flipdim(in,1),2);

% M=2*(N-1)+1;
% out=zeros(M);
% out(N(1):M(1),N(2):M(2))=in;
% out(1:N(1)-1,N(2)+1:M(2))=flipdim(in(2:N(1),2:N(2)),1);
% out(N(1)+1:M(1),1:N(2)-1)=flipdim(in(2:N(1),2:N(2)),2);
% out(1:N(1),1:N(2))=flipdim(flipdim(in(1:N(1),1:N(2)),1),2);
end