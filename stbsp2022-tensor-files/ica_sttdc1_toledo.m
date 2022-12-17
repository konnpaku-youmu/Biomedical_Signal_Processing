function [Fe,FNORMe]=ica_sttdc1(Y,nrs)
% ICA BY SIMULTANEOUS THIRD-ORDER TENSOR DIAGONALIZATION (COMPLEX CASE)
% (VERSION 1)
% 
% [Fe,FNORMe]=ica_sttdc1(Y) performs "Independent Component Analysis" (ICA)
% or "Blind Source Separation", in the following sense. 
%
% The *data model* is given by Y = F X + N, in which: 
%   1. Y is an (nxT)-matrix of observations,
%   2. F is an (nxnrs)-dimensional "transfer" or "mixing" matrix with linearly
%      independent columns,
%   3. X is an (nrsxT)-matrix containing T snapshots of nrs statistically 
%      independent source signals. These signals are assumed to be complex,
%      zero-mean, and at most one of them has zero fourth-order cumulant 
%      (according to the definition Cum_{ijkl} = cum{x_i,x_j^*,x_k^*,x_l}).
%   4. N is an (nxT)-matrix representing additive noise. 
%
% The *goal* is to estimate F and X, when only Y is observed. Denoting the 
% estimates by Fe and Xe, Y is decomposed as Y = Fe Xe.
%
% The *means* on which the algorithm is based, are:
%   1. Prewhitening and projection ("standardization"): the covariance of Y
%      is diagonalized, the dimensionality of the data is reduced to nrs,
%      N is treated as spatially white noise.
%   2. Simultaneous diagonalization of the third-order tensors, obtained by 
%      fixing the fourth index in the fourth-order cumulant Cum_{ijkl} of the 
%      standardized data. N is treated as Gaussian noise.
% Ref: De Lathauwer L., De Moor B., Vandewalle J., ``Independent component analysis 
%		and (simultaneous) third-order tensor diagonalization'', 
%		IEEE Transactions on Signal Processing, vol. 49, no. 10, Oct. 2001, pp. 2262-2271.	
%
% Input:  Y: (nxT)-dataset of observations. When Y is given in (Txn)-format 
%	    the datamodel Y = X*F' + N' is assumed. 
%	  nrs: optional argument that fixes the "number-of-sources". If omitted,
%	       nrs is assumed to be equal to n.
% Output: Fe: "estimate-of-F". The column lengths are such that every estimated
%	      source signal has unit variance.
%	  FNORMe: normalized representation of Fe. Unit-norm columns, with
%                 largest-modulus element real and positif. The ordering 
%		  corresponds to decreasing power of the source estimates.
%
% Copyright
% Noncommercial use only
% Author: Lieven De Lathauwer
% Contact: Lieven.DeLathauwer@kuleuven.be
%

% ( Remarks: this version evaluates the significance of an inner Jacobi-rotation 
% by comparing the modulus of the off-diagonal element to the statistical
% treshold 1/(sqrt(T)*100)) ) 

[nn,TT]=size(Y);T=max(nn,TT);n=min(nn,TT);
if TT==n, Y=Y';end; 			% Y is now (nxT) with n<=T.


%%%% STANDARDIZATION: WHITENING AND PROJECTION

if nargin==1, nrs = n; end;
[U,S,V]=svd(Y',0); 
r=nrs;
if r<n
varn = mean(diag(S(r+1:n,r+1:n)).^2);
else
varn = 0;
end
U=U(:,1:r);
S= diag( sqrt(diag(S(1:r,1:r)).^2 - varn*ones(r,1)) );
V=V(:,1:r);
UFe = V; SFe = S'/sqrt(T); 
Fe = UFe * SFe; Z=inv(SFe)*UFe'*Y;  %%% We have Y = Fe * Z with Fe = UFe*SFe
				    %%% and Z normalized (white; dimension=rxT)

%%%% INITIAL STATISTICS AND CONTRAST

for ii=1:r
  for jj=ii:r
    Cov(ii,jj) = Z(ii,:) * Z(jj,:)'/T;
    Cov(jj,ii) = conj(Cov(ii,jj));
  end
end

for ii=1:r
  for jj=ii:r
    Cov2(ii,jj) = sum(Z(ii,:) .* Z(jj,:)) /T;
    Cov2(jj,ii) = Cov2(ii,jj);
  end
end

for ii=1:r
  for jj=1:r
    for kk=1:r
      for ll=1:r
	M((ii-1)*r+jj,(kk-1)*r+ll) = ...
	        sum(Z(ii,:).*conj(Z(jj,:)).*conj(Z(kk,:)).*Z(ll,:))/T;
	Cum((ii-1)*r+jj,(kk-1)*r+ll) = M((ii-1)*r+jj,(kk-1)*r+ll) ...
		 - Cov(ii,jj)*conj(Cov(kk,ll)) - Cov(ii,kk)*conj(Cov(jj,ll)) ...
		 - Cov2(ii,ll) * conj(Cov2(jj,kk));
      end
    end
  end
end

for ll=1:r
   for kk=1:r
     CUM(:,(ll-1)*r+kk) = Cum(:,(kk-1)*r+ll);
   end
end
% The third-order tensors to be processed are the matrix slices of CUM,
% for which l is fixed. 

%contrast=0; 
%for lt=1:r
%	for jt=1:r
%	contrast = contrast + abs(CUM((jt-1)*r+jt,(lt-1)*r+jt))^2; 
%	end
%end
%fprintf('initial contrast=%g\n',contrast);

%%%% UNITARY TRANSFORMATION BY SIMULTANEOUS THIRD-ORDER TENSOR DIAGONALISATION
 
Rot=eye(r);				      % Rot = overall rotation --> VFe'
nsweep=0; 				      % number of computed sweeps
cntinue = 1; treshold = 1 / (sqrt(T)*100);   % a statistical threshold

while cntinue, cntinue=0;                          %%%%%% start iteration

Q=eye(r);
  for it=1:r-1,
  for j= it+1:r,
    
    M = zeros(3,3);
    
    for ll=1:r
      
      p111 = CUM((it-1)*r+it,(ll-1)*r+it); p112 = CUM((it-1)*r+it,(ll-1)*r+j);
      p121 = CUM((it-1)*r+j,(ll-1)*r+it);  p122 = CUM((it-1)*r+j,(ll-1)*r+j);
      p211 = CUM((j-1)*r+it,(ll-1)*r+it);  p212 = CUM((j-1)*r+it,(ll-1)*r+j);
      p221 = CUM((j-1)*r+j,(ll-1)*r+it);   p222 = CUM((j-1)*r+j,(ll-1)*r+j);
    
      a111 = abs(p111)^2 + abs(p222)^2;
      a112 = p112 + p121;
      a122 = p212 + p221;
      b1 = 0.25 * ( a111 + abs(a112)^2 + abs(p211)^2 + abs(p122)^2 + ... 		
	   abs(a122)^2 ) + 0.5 * real( p111*conj(a122) + p222*conj(a112) );
      b2 = 0.5 * real( p111*conj(p122) + p211*conj(a112) + ...
      		       p222*conj(p211) + p122*conj(a122) );
      
      
      m11 = a111;
      m12 = 0.5 * imag( p111*conj(a112) + p211*conj(p111) + ...
		conj(p222)*p122 + p222*conj(a122) );
      m22 = b1 - b2;
      m31 = 0.5 * real( p222*conj(a122+p122) - p111*conj(a112+p211) );
      m32 = -0.5 * imag( p111*conj(p122) + p211*conj(a112) + ...
		   p211*conj(p222) + conj(p122)*(a122) );
      m33 = b1 + b2;
      m21 = m12;
      m13 = m31;
      m23 = m32;
      
      M = M + [m11 m12 m13; m21 m22 m23; m31 m32 m33];
       
    end
    
    [e,d] = eig(M);
    [maximum,pos] = max(diag(d));
    X = e(:,pos);
    X = X / norm(X);
    if X(1) < 0, X = -X; end	%%% X is now [cos(2t) sin(2t)sin(phi) 			
				%%%			sin(2t)cos(phi)] 
    				%%% with t an inner rotation angle
    				%%% and phi the angle of the complex unit
    
    co = sqrt((X(1)+1)/2);
    si_tot = (X(3)-i*X(2))/(2*co);
    
    apply = abs(si_tot) > treshold;
    cntinue = cntinue | apply;	% check for one significant rotation
    
    if apply
    
    Qij=eye(r); 
    Qij(it,it)=co; Qij(it,j)=-conj(si_tot); Qij(j,it)=si_tot; Qij(j,j)=co;
    
    for lt=1:r
     CUM(:,(lt-1)*r+1:lt*r) = kron(Qij,conj(Qij))*CUM(:,(lt-1)*r+1:lt*r) * Qij'; 
    end
    
    %contrast1=0; 
	%for l2=1:r
	%for j2=1:r
	%contrast1 = contrast1 + abs(CUM((j2-1)*r+j2,(l2-1)*r+j2))^2; 
	%end
	%end
    %fprintf('intermediate contrast1=%g\n',contrast1);
	
    Q=Qij*Q;
    
    end 	% end if
  end;		% end loop it
  end;		% end loop j

nsweep = nsweep+1;  
Rot=Rot*Q';

end;                           %%%%%% stop iteration (end while)

Fe=Fe*Rot;
VFe = Rot';

%%%%%% FINAL CONTRAST
%
%contrast=0; 
%for lt=1:r
%	for jt=1:r
%	contrast = contrast + abs(CUM((jt-1)*r+jt,(lt-1)*r+jt))^2; 
%	end
%end
%fprintf('final contrast =%g\n',contrast);

%%%% NORMALIZATION TO FNORMe

delta=diag(sqrt(sum(Fe.*conj(Fe)))); 	% column norms

[d,I]=sort(-diag(delta));E=eye(r);
P=E(:,I)';delta=P*delta*P';FNORMe=Fe*P'; % descending order

FNORMe=Fe*inv(delta);			% length normalisation

[y,I]=max(abs(FNORMe));
for it=1:r,Lambda(it)=conj(FNORMe(I(it),it));end;Lambda=Lambda./abs(Lambda);
FNORMe=FNORMe*diag(Lambda);		% phase normalisation
