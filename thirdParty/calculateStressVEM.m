function [stress,strain]=calculateStressVEM(G,uu, op,varargin)
% This is a third party SINTEF function slightly modified. 

opt.do_patch=false;
opt=merge_options(opt,varargin{:});

if(G.griddim==2)
    lindim=3;
else
    lindim=6;
end
stress=reshape(op.D*op.WC'*op.assemb'*reshape(uu',[],1),lindim,[])';
if(opt.do_patch)
  stress=patchRecovery(G,stress);
end
if(G.griddim==3)
 stress=bsxfun(@rdivide,stress(:,[1:3,5,6,4]),[ones(1,3),2*ones(1,3)]);
else
    assert(G.griddim==2)
  stress=bsxfun(@rdivide,stress,[ones(1,2),2]);  
end

if(nargout==2)
    %warning('Strains need to be checked')
    strain=reshape(op.WC'*op.assemb'*reshape(uu',[],1),lindim,[])';
end

end