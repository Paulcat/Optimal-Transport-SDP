function S = genorder(n,ordering,positive)
%GENORDER Generate colexicographical ordering
%   C = GENCOLEX(N) produces a sequence of (multi-)indices, starting from
%   (0,...,0) and up to (N1,...,Nd), sorted in the colexicographic order.
%
%   C = GENCOLEX(N,'sym') produces indices ranging from (-N1,...,-Nd) to
%   (N1,...,Nd).


d = numel(n);

switch d
    case 1
        S = (-n:n)';
    case 2
        [Y,X] = meshgrid(-n(1):n(1),-n(2):n(2));
        switch ordering
            case 'colex'
                S = [X(:),Y(:)];
            case 'lex'
                S = [Y(:),X(:)];
            case 'gcolex'
                S = [X(:),Y(:)];
                [~,i] = sort(sum(abs(S),2));
                S = S(i,:);
            case 'glex'
                S = [Y(:),X(:)];
                [~,i] = sort(sum(abs(S),2));
                S = S(i,:);
        end
end

if positive
    S(~all(S>=0,2),:) = [];
end
            


end

