function [A,loc] = integralEbd(varargin)

if nargin == 7
    % --> \int_{mesh(x)} \int_{mesh(y)} psi(x)' f(x,y) psi(y) dxdy
    X    = varargin{1};
    Ydom    = varargin{2};
    green   = varargin{3};
    k       = varargin{4};
    v       = varargin{5};
    tol     = varargin{6};
    lambda  = varargin{7};
    
    
    
    Nx     = size(X,1);
    
    % Domain with quadrature
    [Y,Wy] = Ydom.qud;
    Ny     = size(Y,1);
    Wy     = spdiags(Wy,0,Ny,Ny);
    
    
    % Finite element matrix with integration
    Mv = v.uqm(Ydom);
    if iscell(Mv)
        Mv{1} = Wy * Mv{1};
        Mv{2} = Wy * Mv{2};
        Mv{3} = Wy * Mv{3};
    else
        Mv = Wy * Mv;
    end
    
    [Gxy,loc]       = MVproduct(green,X(:,1:2),Y(:,1:2),tol,k,'lambda',lambda);
    aP              = @(V) Gxy(Mv*V);
    loc             = loc*Mv;
    A               = AbstractMatrix([],aP,Nx,size(Mv,2));
    
    
    %%% EFFCIENT BESSEL DECOMPOSITION WITH BOUNDARY ELEMENT OPERATOR
elseif nargin == 8
    % --> \int_{mesh(x)} \int_{mesh(y)} psi(x)' f(x,y) psi(y) dxdy
    
    Xdom    = varargin{1};
    Ydom    = varargin{2};
    u       = varargin{3};
    green   = varargin{4};
    k       = varargin{5};
    v       = varargin{6};
    tol     = varargin{7};
    lambda  = varargin{8};
    
    
    
    % Domain with quadrature
    [X,Wx] = Xdom.qud;
    
    Nx     = size(X,1);
    Wx     = spdiags(Wx,0,Nx,Nx);
    
    % Domain with quadrature
    [Y,Wy] = Ydom.qud;
    Ny     = size(Y,1);
    Wy     = spdiags(Wy,0,Ny,Ny);
    
    % Finite element matrix with integration
    Mu = u.uqm(Xdom);
    if iscell(Mu)
        Mu{1} = Mu{1}' * Wx;
        Mu{2} = Mu{2}' * Wx;
        Mu{3} = Mu{3}' * Wx;
    else
        Mu = Mu' * Wx;
    end
    
    % Finite element matrix with integration
    Mv = v.uqm(Ydom);
    if iscell(Mv)
        Mv{1} = Wy * Mv{1};
        Mv{2} = Wy * Mv{2};
        Mv{3} = Wy * Mv{3};
    else
        Mv = Wy * Mv;
    end
    
    [Gxy,loc] = MVproduct(green,X(:,1:2),Y(:,1:2),tol,k,'lambda',lambda);
    aP        = @(V) Mu*Gxy(Mv*V);
    loc       = Mu*loc*Mv;
    A         = AbstractMatrix([],aP,size(Mu,1),size(Mv,2));
end


end




