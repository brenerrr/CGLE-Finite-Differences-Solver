function [E] = lhs(A, alpha, theta, p, l, beta)

f = - complex(1, alpha) * theta * p;

nA = length(A);
E = zeros(nA, nA);
for i = 2:nA-1 
    e = 1 + 2 * complex(1, alpha) * theta * p + l * complex(1, beta) * abs(A(i))^2;
    E(i,i) = e;
    E(i,i-1) = f; 
    E(i,i+1) = f;
end

% Boundary conditions 
e1 = 1 + 4/3 * complex(1, alpha) * theta * p + l * complex(1, beta) * abs(A(1))^2;
en = 1 + 4/3 * complex(1, alpha) * theta * p + l * complex(1, beta) * abs(A(end))^2;

E(1,1) = e1;
E(1,2) = 5/6*f;
E(1,end) = 2/3*f;
E(1,end-1) = -1/6*f;

E(end,1) = 2/3*f;
E(end,2) = -1/6*f;
E(end,end) = en;
E(end,end-1) = 5/6*f;
end