function energy = ljp(distance, param)

epsilon = param(1);
sigma = param(2);
r = sigma / distance;
energy = 4*epsilon*(r^12 - r^6);

end

