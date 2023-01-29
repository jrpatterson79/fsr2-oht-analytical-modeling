function [h] = compass3(col_vecs)

num_cols = size(col_vecs,2);
bnds = [ min(col_vecs'); max(col_vecs')];

figure
hold on
for i = 1:1:num_cols
    plot3([0 col_vecs(1,i)],[0 col_vecs(2,i)],[0 col_vecs(3,i)],'b');
end
plot3([0; 0],[0; 0],[0; 1],'-.r');
plot3([0; 0],[0; 1],[0; 0],'-.r');
plot3([0; 1],[0; 0],[0; 0],'-.r');

hold off
axis vis3d
axis equal
grid on
