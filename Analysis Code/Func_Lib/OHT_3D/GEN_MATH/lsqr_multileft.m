function xsol = lsqr_multileft(A,b)

setup.type = 'crout';
setup.milu = 'row';
setup.droptol = 0.05;
[L,U] = ilu(A,setup);
num_left = size(b,2);
num_unk = size(b,1);
xsol = zeros(num_unk,num_left);
parfor i = 1:1:num_left
    xsol(:,i) = lsqr(A,b(:,i),1e-4,600,L,U,[]);
end