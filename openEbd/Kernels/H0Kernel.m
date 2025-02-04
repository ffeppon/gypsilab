classdef H0Kernel < Kernel
    properties
        k; 
        C; % Such that G(r) = C*besselh(0,1,k*r)
        % Hankel function of the first kind. 
        cutoff
        val0
    end
    methods
        function[this] = H0Kernel(kk,CC,cut,valz)
            if ~exist('cut','var')||isempty(cut)
                cut = 0;
            end
            if ~exist('valz','var')||isempty(valz)
                valz = 0;
            end
            if ~exist('CC','var')||isempty(CC)
                CC = 1;
            end
            if ~exist('kk','var')||isempty(kk)
                kk = 1;
            end
            this.cutoff = cut; 
            this.val0 = valz;
            this.k = kk;
            this.C = CC;
            Jk = J0Kernel(kk);
            Yk = Y0Kernel(kk);
            modelKern = this.C*(Jk + 1i*Yk);
            this.func = @(x)(H0Kernel.funcCutoff(this.k,x,this.cutoff,this.val0));
            this.der = modelKern.der;
            this.scalFunc = modelKern.scalFunc;
        end
        function[out] = dilatation(this,lambda)
            cut = this.cutoff/lambda;
            valz = this.val0;
            out = H0Kernel(this.k*lambda,this.C,cut,valz);
        end
        function[out] = mtimes(this,mu)
            if isa(this,'Kernel')
                assert(and(isa(mu,'double'),isscalar(mu)));
                out = H0Kernel(this.k,this.C*mu);
            else
                out = mtimes(mu,this);
            end
        end
        function[rq] = radialQuadKernel(this,a,tol,varargin)
            modelKern = J0Kernel(this.k) + 1i*Y0Kernel(this.k);
            rqTemp = modelKern.radialQuadKernel(a,tol/abs(this.C),varargin{:});
            rq = this.C*rqTemp;
        end
    end
    methods (Static)
        function[y] = funcCutoff(k,x,cut,valz)
            y = besselh(0,k*x);
            y(x < cut) = valz;
        end
    end
end
