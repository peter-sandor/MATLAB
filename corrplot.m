function corrplot(x)
a=load(['trace' num2str(x)]);
subplot(211);
plot(a);
% b=imread(['pic' num2str(x) '.jpg']);
b=monopic(imread(['pic' num2str(x) '.jpg']));
subplot(212);
imagesc(b);
colorbar
end