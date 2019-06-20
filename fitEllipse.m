function [parameters, coefficients] = fitEllipse(X, Y)

% normalize data
mx = mean(X);
my = mean(Y);
sx = (max(X)-min(X))/2;
sy = (max(Y)-min(Y))/2; 

x = (X-mx)/sx;
y = (Y-my)/sy;

% Force to column vectors
x = x(:);
y = y(:);

% Build design matrix
D = [ x.*x  x.*y  y.*y  x  y  ones(size(x)) ];

% Build scatter matrix
S = D'*D;

% Build 6x6 constraint matrix
C(6,6) = 0; C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;

% Solve eigensystem
[gevec, geval] = eig(S,C);

% Find the negative eigenvalue
neg_eigval = find(real(diag(geval)) < 1e-8 & ~isinf(diag(geval)));

% Extract eigenvector corresponding to negative eigenvalue
A = real(gevec(:,neg_eigval));

% unnormalize
coefficients = [
  A(1)*sy*sy,   ...
      A(2)*sx*sy,   ...
      A(3)*sx*sx,   ...
      -2*A(1)*sy*sy*mx - A(2)*sx*sy*my + A(4)*sx*sy*sy,   ...
      -A(2)*sx*sy*mx - 2*A(3)*sx*sx*my + A(5)*sx*sx*sy,   ...
      A(1)*sy*sy*mx*mx + A(2)*sx*sy*mx*my + A(3)*sx*sx*my*my   ...
      - A(4)*sx*sy*sy*mx - A(5)*sx*sx*sy*my   ...
      + A(6)*sx*sx*sy*sy   ...
      ]';

%% Convert conic coefficients to ellipse parameters
% Convert to geometric radii, and centers
thetarad = 0.5*atan2(coefficients(2),coefficients(1) - coefficients(3));

cost = cos(thetarad);
sint = sin(thetarad);
sin_squared = sint.*sint;
cos_squared = cost.*cost;
cos_sin = sint .* cost;

Ao = coefficients(6);
Au =   coefficients(4) .* cost + coefficients(5) .* sint;
Av = - coefficients(4) .* sint + coefficients(5) .* cost;
Auu = coefficients(1) .* cos_squared + coefficients(3) .* sin_squared + coefficients(2) .* cos_sin;
Avv = coefficients(1) .* sin_squared + coefficients(3) .* cos_squared - coefficients(2) .* cos_sin;

tuCentre = - Au./(2.*Auu);
tvCentre = - Av./(2.*Avv);
wCentre = Ao - Auu.*tuCentre.*tuCentre - Avv.*tvCentre.*tvCentre;

uCentre = tuCentre .* cost - tvCentre .* sint;
vCentre = tuCentre .* sint + tvCentre .* cost;

Ru = -wCentre./Auu;
Rv = -wCentre./Avv;

Ru = sqrt(abs(Ru)).*sign(Ru);
Rv = sqrt(abs(Rv)).*sign(Rv);

parameters = [uCentre, vCentre, Ru, Rv, thetarad];

end

