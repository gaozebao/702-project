function [V,d] = getbasis(B)
  [n,m] = size(B);
  [V,d] = helper(B);
  [vn,p] = size(V);
  assert(vn==n);
  V = [ones(n,1), V];
  [V,i] = reorderbasis(V);
  d = d(i);
  
  i = find(d~=0);
  d = d(i);
  V = V(:,i);
  
  B1 = V * diag(d) * V';
  err = sum(sum(abs(B-B1)));
  if err > 1e-6
      V,d,B1,B
  end
  assert(err <= 1e-6);
  
function [V,d] = helper(B)
  tol = 1e-8;
  [n,m] = size(B);
  d0 = min(min(B));
  B = B - d0 .* ones(n,n);
  [i,j] = find(B == 0, 1, 'first');
  left = (B(i,:) > tol)';
  right = ones(n,1) - left;
  [leftV,leftd] = subtreehelper(B,left);
  [rightV,rightd] = subtreehelper(B,right);
  
  V = [left leftV right rightV];
  d = [d0;leftd;rightd];
  
function [V,d] = subtreehelper(B,v)
  [n,m] = size(B);
  seq = 1:n;
  nv = sum(v==1);
  if (nv == 1)
    V = [];
    d = diag(B);
    d = d(v==1);
    return;
  end
  
  ind = find(v);
  [subV,d] = helper(B(ind,ind));
  [vn,p] = size(subV);
  V = zeros(n,p);
  V(ind,:) = subV;
  
  
