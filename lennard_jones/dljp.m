function denergy = dljp(distance, param)

epsilon = param(1);
sigma = param(2);
r = sigma / distance;
% energy = 4*epsilon*(r^12 - r^6);
denergy = 4 * epsilon * (-12 / distance * r^12 + 6 / distance * r^6);

end