function corrplot(x)
a=load(['trace' num2str(x)]);
subplot(211);
plot(a);
b=imread(['pic' num2str(x) '_0.jpg']);
subplot(212);
imagesc(b);
colorbar
end