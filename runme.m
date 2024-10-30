%-------------------------------------------------------
% example script for compute riemann mapping from a triangulated polygonal
% domain to a triangle. Input is a planar mesh and 3 boundary points given
% in counter clockwise orientation
%
% Code written by Nadav Dym and Raz Slutsky, based on previous code by Noam Eigerman. inquiries about the code can be sent to
% nadavdym@gmail.com
% Please cite
% "A Linear Variational Principle for Riemann Mappings and Discrete
% Conformality" Dym, Lipman, Slutsky
% and
% "Orbifold Tutte embeddings" Aigerman and Lipman
%-------------------------------------------------------

load('triangulated_koch.mat');
[V_flat,data]=discrete_conformal(V,T,chosen_points);
visualize_mapping(V,V_flat,T,data);
