function [V_flat,data]=discrete_conformal(V,F,Q)
%----------------------------------------------------------------------------------------------------
% compute discrete conformal mapping from mesh defined by V,F with three
% selected points Q to a isosceles triangle
%
% Code written by Nadav Dym and Raz Slutsky. inquiries about the code can be sent to nadavdym@gmail.com
% Please cite
% "A Linear Variational Principle for Riemann Mappings and Discrete Conformality" Dym, Lipman, Slutsky
% and
% "Orbifold Tutte embeddings" Aigerman and Lipman
%-----------------------------------------------------------------------------------------------------



n_v=max(size(V));
P=[0 0;1 0;0 1]'; %vertices of the triangle we map to
%-----------------------------------------------------------------------
% compute boundary and 3point inds
%-----------------------------------------------------------------------
%set boundary conditions
V(:,3)=zeros(size(V,1),1);
tri=triangulation(F,V(:,1),V(:,2));
boundary=tri.freeBoundary();
boundary=boundary(:,1);

%select boundary points
three_point_inds=nan(1,3);
for i = 1:3
    three_point_inds(i)=knnsearch(V(boundary,:),Q(i,:));
end

check_counterclockwise(three_point_inds);


n_boundary=length(boundary);
n_con=n_boundary+3; %cones have double equality constraints
%-----------------------------------------------------------------------
% compute laplacian
%-----------------------------------------------------------------------
L=cotmatrix(V,F);

%check if delaunay
Lnodiag=L-diag(diag(L));
if min(Lnodiag(:))<0
    fprintf('mesh is not delaunay');
end

%-----------------------------------------------------------------------
%get linear constraints
%-----------------------------------------------------------------------
E=[1 2;2 3;3 1]';
b=sparse(2*n_v+n_con,1);
all_x_ind=[];
all_y_ind=[];
val=[];
current_const_num=0;
for ii=1:3
    [c,b_loc]=get_line(P(:,E(1,ii)),P(:,E(2,ii)));
    ind=get_ind(n_boundary,three_point_inds(E(1,ii)),three_point_inds(E(2,ii)));
    new_const_num=length(ind);
    new_const_ind=current_const_num+1:current_const_num+new_const_num;
    all_x_ind=[all_x_ind; new_const_ind'; new_const_ind'];
    all_y_ind=[all_y_ind; boundary(ind); n_v+boundary(ind)];
    val=[val; c(1)*ones(size(ind)); c(2)*ones(size(ind))];
    b(2*n_v+new_const_ind)=b_loc;
    current_const_num=current_const_num+new_const_num;
end
A=sparse(all_x_ind,all_y_ind,val,n_con,2*n_v);
%-----------------------------------------------------------------------
% solve linear system
%-----------------------------------------------------------------------
H=[blkdiag(L,L) A';
    A sparse(n_con,n_con)];
t=tic;
x=H\b;
lin_time=toc(t);
fprintf('linsolver time %d',lin_time);
x=full(x);
V_flat=[x(1:n_v) x(n_v+1:2*n_v) zeros(n_v,1)];

data.lintime=lin_time;
data.boundary=boundary;
data.three_point_inds=three_point_inds;

end


function [c,b]=get_line(p,q)
%get line between the points p and q
v=p-q;
n=[v(2); -v(1)];
c=n/norm(n);
b=c'*p;
end

function ind=get_ind(n,first,last)
%get boundary indices between points first and last
if first<last
    ind=first:last;
else
    ind=[first:n 1:last];
end
ind=ind';
end

function check_counterclockwise(three_point_inds)
%check input points are orientated correctly
temp=[three_point_inds three_point_inds];
[~,min_ind]=min(three_point_inds);
final_inds=temp(min_ind:min_ind+2);
if norm(final_inds-sort(final_inds,'ascend'))>0
    error('input three points not oriented counterclockwise')
end
end













