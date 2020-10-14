function F = TDSE_movie(solution_in)
fig1=figure('Visible','off');
ax1=subplot(211);
ax2=subplot(212);
fprintf(1, 'Creating movie:\n');
for ind1=1:length(solution_in.t)
    plot(ax1,solution_in.x,[abs(solution_in.Psi(:,ind1)) real(solution_in.Psi(:,ind1))]);
    xlim(ax1,[min(solution_in.x) max(solution_in.x)]);
    ylim(ax1,[-max(max(real(solution_in.Psi))) max(max(real(solution_in.Psi)))]);
    legend(ax1,'|\Psi|','Re(\Psi)')
    title(ax1,['time = ' num2str(roundP(solution_in.t(ind1),2)) ' [at.u.]'])
    plot(ax2,solution_in.x,real(solution_in.V(:,ind1)),'k-');
	xlim(ax2,[min(solution_in.x) max(solution_in.x)]);
    if min(min(solution_in.V)) ~= max(max(solution_in.V))
        ylim(ax2,[min(min(real(solution_in.V))) max(max(real(solution_in.V)))]);
    end
    xlabel(ax2,'distance [at.u.]');
    legend(ax2,'Potential')
    F(ind1) = getframe(fig1);
    fprintf(1, '.');
    if mod(ind1,100)==0
        fprintf(1, '\n');
    end
end
fprintf(1, '\n');
disp('Done.');
end