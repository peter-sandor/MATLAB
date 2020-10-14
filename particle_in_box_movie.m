% Movie

N_expand = 29;
n = (1:N_expand)';
a = 1/(sqrt(2)*pi).*(1./(n./2-1).*sin((n./2-1)*pi)-1./(n./2+1).*sin((n./2+1)*pi));

a(2) = 1/sqrt(2);
for k = 4:2:N_expand-1
a(k)= 0;
end

h = 6.636e-34; % [m^2*kg/s)]
mass = 9.11e-31; % [kg]
L = 1e-10; % [m]
N_frame = 100;
x = 0:0.01:1;
omega = h/(2*pi)*pi^2*n.^2/(8*mass*L);
T = 2*pi/omega(1);
amp = zeros([length(x),length(n)]);
h_fig = figure;
for j = 1:N_frame
    t = T*j/N_frame;
    for k = 1:length(n)
        amp(:,k) = a(k)*sin(n(k)*pi*x)*real(exp(-i*omega(k)*t));
    end
    psi = zeros([length(x),1]);
    psi = (sum(amp(:,1:N_expand),2)).^2;
    
    plot(x,psi,'k.')
%     hold on; plot(x,amp(:,1),'y')
%     plot(x,amp(:,2), 'r')
%     plot(x,amp(:,3), 'b')
%     plot(x,amp(:,5), 'g')
%     plot(x,amp(:,7), 'm')
%     hold off;
    F(j) = getframe(h_fig);
end
close(h_fig)
% h_fig2 = figure;
% movie(h_fig2,F,10,25)
movie2avi(F,'box3.avi','compression','Cinepak');