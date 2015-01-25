function [ q1 ] = projectQperp( Q,q )

%   Detailed explanation goes here
%   q1 is normalized
for i=1:size(Q,2)
    q=q-Q(:,i)*Q(:,i)'*q;
end
q1=q/norm(q);
end

