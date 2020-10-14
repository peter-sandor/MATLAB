Nlines=20;

circle_obj(1).center=[0 0];% [x0 y0]
circle_obj(1).R=sqrt(2);

% circle_obj(2).center=[0.8 0];% [x0 y0]
% circle_obj(2).R=0.3;

line_obj0(1).A=0; % [A B] for y=A*x+B; if A=inf, then line_str(2)=x coord.
line_obj0(1).B=-0.1;
line_obj0(1).k=[cos(atan(line_obj0(1).A)); sin(atan(line_obj0(1).A))]; % unit vector showing direction
line_obj0(1).pstart=[-sqrt(2); line_obj0(1).B]; % start point
line_obj0(1).pend=[sqrt(2); line_obj0(1).B]; % start point

% line_obj0(2).A=-3; % [A B] for y=A*x+B; if A=inf, then line_str(2)=x coord.
% line_obj0(2).B=-1;
% line_obj0(2).k=[cos(atan(line_obj0(2).A)); sin(atan(line_obj0(2).A))]; % unit vector showing direction
% line_obj0(2).pstart=[-sqrt(2); line_obj0(1).B]; % start point
% line_obj0(2).pend=[sqrt(2); line_obj0(1).B]; % start point

line_obj(1).A=0; % [A B] for y=A*x+B; if A=inf, then line_str(2)=x coord.
line_obj(1).B=0.1;
line_obj(1).k=[cos(atan(line_obj(1).A)); sin(atan(line_obj(1).A))]; % unit vector showing direction
line_obj(1).pstart=[0 -0.05]; % start point
line_obj(1).pend=[]; % end point
%%
ind1=1;
while ind1<=Nlines
    [object1,point1]=next_object(line_obj(ind1),line_obj0,circle_obj);
    line_obj(ind1).pend=point1;
    line_obj(ind1+1)=calc_reflection(line_obj(ind1),object1,point1);
%     disp([num2str(ind1) '/' num2str(Nlines)]);
    ind1=ind1+1;
end
%%
Ndraw=100;
figure;axes;hold on;

for ind1=1:length(circle_obj)
    xvec=circle_obj(ind1).center(1)-circle_obj(ind1).R:2*circle_obj(ind1).R/(Ndraw-1):circle_obj(ind1).center(1)+circle_obj(ind1).R;
    line(xvec,sqrt(circle_obj(ind1).R^2-(xvec-circle_obj(ind1).center(1)).^2)+circle_obj(ind1).center(2),'color','k');
    line(xvec,-sqrt(circle_obj(ind1).R^2-(xvec-circle_obj(ind1).center(1)).^2)+circle_obj(ind1).center(2),'color','k');
end

for ind1=1:length(line_obj0)
    xstep=(line_obj0(ind1).pend(1)-line_obj0(ind1).pstart(1))/(Ndraw-1);
    xvec=map2colvec(line_obj0(ind1).pstart(1):xstep:line_obj0(ind1).pend(1));
    lines0(:,:,ind1)=[xvec line_obj0(ind1).A*xvec+line_obj0(ind1).B];
    line(lines0(:,1,ind1),lines0(:,2,ind1),'color','k','linewidth',1);
end

for ind1=1:length(line_obj)
    xstep=(line_obj(ind1).pend(1)-line_obj(ind1).pstart(1))/(Ndraw-1);
    xvec=map2colvec(line_obj(ind1).pstart(1):xstep:line_obj(ind1).pend(1));
    lines1(:,:,ind1)=[xvec line_obj(ind1).A*xvec+line_obj(ind1).B];
    line(lines1(:,1,ind1),lines1(:,2,ind1),'color','r','linewidth',2);
end

xlim([-circle_obj(1).R circle_obj(1).R]);
ylim([-circle_obj(1).R circle_obj(1).R]);
axis square;