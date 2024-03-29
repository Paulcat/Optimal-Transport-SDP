function [J] = Jaccard(x0,xe,r)
%JACCARD_INDEX Computes Jaccard index between two vectors x0 and xe
%   We follow the code described in the thesis of Q. Denoyelle, in section
%   4.3.2
%
%   r: tolerance on the accuracy of detection

detected_spikes = [];

while size(x0,1) > 0 && size(xe,1) > 0
	% Compute, sort and store distances d(x0_i, xe_j), for all i, j
	% compute all indices combinations
	[i,j] = meshgrid( 1:size(x0,1), 1:size(xe,1) );

	% compute all norms ||x0_i - xe_j||
	norms = sqrt(sum( (x0(i(:),:) - xe(j(:),:)).^2, 2 ));
	[test, id] = min(norms);

	if test > r
		break;
	end

	detected_spikes = [detected_spikes; x0(i(id),:)];
	x0(i(id),:)     = [];
	xe(j(id),:)     = [];
end

TP = size(detected_spikes,1);
FP = size(xe,1);
FN = size(x0,1);

J = TP / (TP+FP+FN);


end
