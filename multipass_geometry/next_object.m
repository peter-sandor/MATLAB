function [object_out point_out] = next_object(line_in,line_objects,circle_objects)

PREC=1e-10;
Nline=length(line_objects);
for ind1=1:Nline
    p1=calc_intersection(line_in,line_objects(ind1));
    if ~isempty(p1)
         % check if intersection is in the direction of the ray in question
         % or distance from starting point is not too small
        if ((p1(1)-line_in.pstart(1))*line_in.k(1)+(p1(2)-line_in.pstart(2))*line_in.k(2))<0 ...
                || sqrt((p1(1)-line_in.pstart(1))^2 + (p1(2)-line_in.pstart(2))^2)<PREC
            points(ind1,:)=[NaN;NaN];
        else
            points(ind1,:)=p1;
        end
    else
        points(ind1,:)=[Inf;Inf];
    end
    dist1(ind1)=sqrt((line_in.pstart(1)-points(ind1,1))^2+(line_in.pstart(2)-points(ind1,2))^2);
end    

Ncircle=length(circle_objects);
for ind1=1:Ncircle
    p2=calc_intersection(line_in,circle_objects(ind1));
    if ~isempty(p2)
        for ind2=1:size(p2,1)
            if ((p2(ind2,1)-line_in.pstart(1))*line_in.k(1)+(p2(ind2,2)-line_in.pstart(2))*line_in.k(2))<0 ...
                || sqrt((p2(ind2,1)-line_in.pstart(1))^2 + (p2(ind2,2)-line_in.pstart(2))^2)<PREC
                temp(ind2,:)=[NaN;NaN];
            else
                temp(ind2,:)=p2(ind2,:);
            end  
%             if ((p2(ind2,1)-line_in.pstart(1))*line_in.k(1)+(p2(ind2,2)-line_in.pstart(2))*line_in.k(2))>0 % intersection is in the direction of the ray in question
%                 temp(ind2,:)=p2(ind2,:);
%             else
%                 temp(ind2,:)=[NaN;NaN];
%             end
        end
        dist2=sqrt((line_in.pstart(1)-temp(:,1)).^2+(line_in.pstart(2)-temp(:,2)).^2);
        [sorted2,index2]=sort(dist2,'ascend');
        points(Nline+ind1,:)=temp(index2(1),:);
    else
        points(Nline+ind1,:)=[Inf;Inf];
    end
    dist1(Nline+ind1)=sqrt((line_in.pstart(1)-points(Nline+ind1,1))^2+(line_in.pstart(2)-points(Nline+ind1,2))^2);
end

% ind_del=vec2ind(dist1<PREC);
% if ~isempty(ind_del)
%     dist1(ind_del)=[];
%     points(ind_del,:)=[];
%     if ind_del<=Nline
%         line_objects(ind_del)=[];
%         Nline=Nline-1;
%     elseif ind_del>Nline
%         circle_objects(ind_del-Nline)=[];
%     end
% end

[sorted1,index1]=sort(dist1,'ascend');
point_out=points(index1(1),:);
if index1(1)<=Nline
    object_out=line_objects(index1(1));
elseif index1(1)>Nline
    object_out=circle_objects(index1(1)-Nline);
else
    object_out=[];
    point_out=[];
end
end