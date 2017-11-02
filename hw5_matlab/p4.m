load topodata;



lam = rbffit(x,y,z);
[ux,uy] = meshgrid( linspace(0,261,100), linspace(0,399,153));
f = rbfval(lam, x, y, ux, uy);
hold on
%surf(ux, uy, f)
%scatter3(x,y,z,40,'.')
%colorbar
%shading interp;
%colormap(autumn);
%lighting phong;
%camlight right;
%xlabel("x-direction");
%ylabel("y-direction");
%zlabel("Height");
%title("2-ravine drainage");
%hold off

hold on
contour(ux, uy, f, 14)
scatter(x,y, '.')
colorbar
xlabel("x-direction");
ylabel("y-direction");
colormap(winter);
hold off