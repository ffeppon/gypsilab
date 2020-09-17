function mesh = mshClean(mesh,dst)

% Unify duplicate vertex, remove unnecessary vertices and relabel elements
% accordingly.

if ~exist('dst','var')||isempty(dst)
    [~,I,J] = unique(single(mesh.vtx),'rows','stable');
else
    % Range search for specified close distance 
    [J,D] = rangeSearch(mesh.vtx,mesh.vtx,dst);
    M     = rangeMatrix(J,D);
    M     = (M>0);
    
    % Security
    if ~issymmetric(M)
        error('mshClean.m : unavailable case');
    end
        
    % Select pairs for descent sort
    M = triu(M);
    for i = 1:size(M,1)
        % Initialization
        J    = [];
        Jnew = 1;
        
        % Iterate on new neighbours
        while (length(Jnew) > length(J))
            % Column indices
            J = find(M(i,:));
            J = J(J>i);
            
            % Row indices
            [I,~] = find(M(:,J));
            I     = unique(I);
            I     = I(I>i);
            
            % Add new values
            M(i,:) = M(i,:) + sum(M(I,:),1);
            M(I,:) = 0;
            
            % Verify if new entries
            Jnew = find(M(i,:));
            Jnew = Jnew(Jnew>i);
        end
    end
    
    % Security
    if (nnz(M) ~= size(M,1))
        error('mshClean.m : unavailable case');
    end
    
    % Find interactions
    [idx,jdx] = find(M);

    % Unicity
    ind       = (1:size(mesh.vtx,1))';
    ind(jdx)  = idx ;
    [~,I,J]   = unique(ind,'stable');
end

% Update mesh
mesh.vtx = mesh.vtx(I,:);
if (size(mesh.elt,1) == 1)
    J = J';
end
mesh.elt = J(mesh.elt);

% Extract vertex table from element
Ivtx              = zeros(size(mesh.vtx,1),1);
Ivtx(mesh.elt(:)) = 1;
mesh.vtx          = mesh.vtx(logical(Ivtx),:);

% Reorder elements
Ivtx(Ivtx==1) = 1:sum(Ivtx,1);
if size(mesh.elt,1) == 1
    mesh.elt = Ivtx(mesh.elt)';
else
    mesh.elt = Ivtx(mesh.elt);
end

% Empty mesh
if (size(mesh.elt,2) == 0)
    I = 0;
    
% Particles mesh
elseif (size(mesh.elt,2) == 1)
    [~,ind] = unique(mesh.elt);
    I       = ones(size(mesh.elt,1),1);
    I(ind)  = 0;

% Edge mesh
elseif (size(mesh.elt,2) == 2)
    I = (mesh.elt(:,1)==mesh.elt(:,2));
    
% Triangular mesh
elseif (size(mesh.elt,2) == 3)
    I = (mesh.elt(:,1)==mesh.elt(:,2)) + ...
        (mesh.elt(:,1)==mesh.elt(:,3)) + ...
        (mesh.elt(:,2)==mesh.elt(:,3)) ;

% Tetrahedral mesh
elseif (size(mesh.elt,2) == 4)
    I = (mesh.elt(:,1)==mesh.elt(:,2)) + ...
        (mesh.elt(:,1)==mesh.elt(:,3)) + ...
        (mesh.elt(:,1)==mesh.elt(:,4)) + ...
        (mesh.elt(:,2)==mesh.elt(:,3)) + ...
        (mesh.elt(:,2)==mesh.elt(:,4)) + ...
        (mesh.elt(:,3)==mesh.elt(:,4)) ;

% Unknown type
else
    error('mshClean.m : unavailable case')
end

% Extract non degenerated elements
if sum(I)
    mesh = mesh.sub(~I);
end
end


function M = rangeMatrix(J,D)
% Row indices and distances
idx = cell2mat(J')';
val = cell2mat(D')'+eps;

% Column indices
jdx = zeros(length(idx)+1,1);
n   = 1;
for i = 1:length(J)
    jdx(n) = jdx(n) + 1;
    n      = n + length(J{i});
end
jdx = cumsum(jdx(1:end-1));

% Sparse Matrix
M = sparse(idx,jdx,val);
end
