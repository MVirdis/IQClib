function [psi,Anu,Bnu,Cnu,Dnu] = basisTF(nu,rho,nq)
%basisTF(nu,rho) makes a vector of fixed basis functions with McMillan degree nu
%and pole location rho

psi = tf(1)*zeros((nu+1),1);

psi(1) = tf(1);

poly_ = [1,-rho]; % base polynomial
poly = poly_;
for iTF=2:nu+1
    psi(iTF,1) = tf(1,poly);
    poly = conv(poly, poly_);
end

% Inflate according to nq (kronecker product)
psi_ = tf(1)*zeros((nu+1)*nq,nq);
psi_(1:nq,1:nq) = tf(1)*eye(nq);
for iTF=1:nu
    psi_((iTF)*nq+1:(iTF+1)*nq,1:nq) = psi(iTF+1)*eye(nq);
end
psi = psi_;

% Extract ss
psi_ss = ss(psi);
Anu = psi_ss.A;
Bnu = psi_ss.B;
Cnu = psi_ss.C;
Dnu = psi_ss.D;

end
