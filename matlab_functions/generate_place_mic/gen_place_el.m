function place_el_rectang = gen_place_el(N_x,M_y,d_x,d_y,Pl)

    place_el_rectang = zeros(N_x*M_y,3);
    if (Pl == 1)
        for i = 1:M_y
            for k = 1:N_x
                place_el_rectang((i-1)*N_x + (k),:) = [0 (-(N_x/2)+k-1)*d_x+d_x/2 ((M_y/2)-i+1)*d_y-d_y/2 ];
            end
        end
    end
    if (Pl == 2)
        for i = 1:M_y
            for k = 1:N_x
                place_el_rectang((i-1)*N_x + (k),:) = [(-(N_x/2)+k-1)*d_x+d_x/2 0 ((M_y/2)-i+1)*d_y-d_y/2 ];
            end
        end
    end
    if (Pl == 3)
        for i = 1:M_y
            for k = 1:N_x
                place_el_rectang((i-1)*N_x + (k),:) = [(-(N_x/2)+k-1)*d_x+d_x/2 ((M_y/2)-i+1)*d_y-d_y/2 0];
            end
        end
    end
    

    
end

