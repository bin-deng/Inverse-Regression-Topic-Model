function p = projectPhi(phi_star, phi, delta)
p = phi_star;
n = length(phi_star);
for i = 1:n
    if (phi(i) == 0 && phi_star(i) < -delta(i))
        p(i) = -delta(i);
    end
    if (phi(i) == 0 && phi_star(i) > delta(i))
        p(i) = delta(i);
    end
    if (phi_star(i) < phi(i)-delta(i))
        if (phi(i) > 0 && phi(i) >= delta(i))
            p(i) = phi(i)-delta(i);
        end
        if (phi(i) < 0 && phi(i)+ delta(i) <= 0 )
            p(i) = phi(i)-delta(i);
        end
        if (phi(i) < 0 && phi(i) > delta(i))
            p(i) = phi(i)-delta(i);
        end
    end
    if (phi_star(i) > phi(i)+delta(i))
        if (phi(i) > 0 && phi(i)-delta(i) >= 0)
            p(i) = phi(i)+delta(i);
        end
         if (phi(i) > 0 && phi(i)-delta(i) < 0)
            p(i) = phi(i)+delta(i);
         end
         if (phi(i) < 0 && phi(i)+delta(i) <= 0)
            p(i) = phi(i)+delta(i);
        end
    end
    if (phi(i) > 0 && phi(i)-delta(i)<0 && phi_star(i) < 0)
        p(i) = 0;
    end
    if (phi(i) < 0 && phi(i)+delta(i) > 0 && phi_star(i) > 0)
        p(i) = 0;
    end
end


end
