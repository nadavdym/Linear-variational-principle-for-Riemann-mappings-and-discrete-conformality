function visualize_mapping(V,V_flat,T,data)
%------------------------------------------------------------------------------------------------------%
% visualize the computed discrete Riemann mappping
% 
% Code written by Nadav Dym and Raz Slutsky. inquiries about the code can be sent to nadavdym@gmail.com
% Please cite
% "A Linear Variational Principle for Riemann Mappings and Discrete Conformality" Dym, Lipman, Slutsky
% and
% "Orbifold Tutte embeddings" Aigerman and Lipman
%------------------------------------------------------------------------------------------------------%
    cone_colors=[1 0.8 0;0.7 0 1; 0 0.5 0.8;0 0 0.5];
    cut_colors=[1 0 0;0 1 0;0 0 1];
    V(:,3)=0;
    b=data.boundary;
    boundary_num=length(b);
    three_point_inds=data.three_point_inds;
    three_point_cyclic=[three_point_inds three_point_inds(1)];
    cones=b(three_point_inds);
    color=V(:,3);
    figure;
    % === Visualization of the original 3D mesh + cones & cuts on it ===
    subplot(1,2,1);
    %draw the mesh
    patch('faces',T,'vertices',V,'facecolor','interp','FaceVertexCData',color);
    hold on;
    %draw the cuts...
    for i=1:3
        path_ind=get_ind(boundary_num,three_point_cyclic(i),three_point_cyclic(i+1));
        curPath=b(path_ind);
        %draw the cut
        line(V(curPath,1),V(curPath,2),'color',cut_colors(i,:),'linewidth',2);
    end
    %draw the cones
    for i=1:length(cones)
        scatter3(V(cones(i),1),V(cones(i),2),V(cones(i),3),100,cone_colors(i,:),'fill');
    end
    
    axis equal
    
    % === Visualization of the flattening + cones & cuts on it ===
    subplot(1,2,2);
    %draw the flattened mesh
    patch('faces',T,'vertices',V_flat,'facecolor','interp','FaceVertexCData',color);
    hold on;
    for i=1:3
        path_ind=get_ind(boundary_num,three_point_cyclic(i),three_point_cyclic(i+1));
        curPath=b(path_ind);
        %draw the cut
        line(V_flat(curPath,1),V_flat(curPath,2),'color',cut_colors(i,:),'linewidth',5);
    end
    %draw the cones
    for i=1:length(cones)
        flat_cones=cones(i);
        scatter(V_flat(flat_cones,1),V_flat(flat_cones,2),100,cone_colors(i,:),'fill');
    end
    axis equal


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