using PlotlyJS
tmp_xr = [
            1 1 1 1 1 1
            1 0 1 1 0 0;
            1 0 1 1 0 1;
            1 0 0 0 0 1;
            1 1 0 0 0 0;
            1 1 1 0 0 0;
            1 1 1 0 0 0;
        ]

tmp_r = [6.0; 5.0; 4; 3;2; 1;0.0]
tmp_x = [5.0;4.0; 3;2; 1;0.0]

function DetermineBoundary(XR,Rv,Xv)
    r = Vector{Float64}()
    x = Vector{Float64}()
    flag_dir = "left"
    flag_end = false
    ind = [(-1,0);(0,1);(1,0);(0,-1)]

    XR = reverse(XR)'
    Rv = reverse(Rv)
    Xv = reverse(Xv)

    lenx = size(XR)[1]
    lenr = size(XR)[2]

    #starting indices
    indr = 1
    indx = findfirst(k->k==1,XR[:,1])

    #saving first entry
    push!(r,Rv[indr])
    push!(x,Xv[indx])

    function next_ind(start_ind,r,x)
        flag_end == false
        for i=start_ind:start_ind+3
            if mod(i,4) == 0
                dx = ind[4][1]
                dr = ind[4][2]
            else
                dx = ind[mod(i,4)][1]
                dr = ind[mod(i,4)][2]
            end

            if indr+dr == 0 || indr+dr > lenr || indx+dx > lenx 
                #skip boundarys
            elseif indx + dx == 1 
                indr = indr + dr
                indx = indx + dx
                push!(r,Rv[indr])
                push!(x,Xv[indx])
                indr = indr -1
                
                while indr != 0 && XR[indx,indr] == 1  # go x-axis to zero 
                    push!(r,Rv[indr])
                    push!(x,Xv[indx])
                    indr = indr -1
                end
                flag_end = true
                break;
            elseif XR[indx+dx,indr+dr] != 1 
                #skip unstable case
            elseif XR[indx + dx,indr + dr] == 1
                indr = indr + dr
                indx = indx + dx
                push!(r,Rv[indr])
                push!(x,Xv[indx])

                if mod(i,4) == 1
                   flag_dir = "bottom"
                   break;
                elseif mod(i,4) == 2
                   flag_dir = "left"
                   break;
                elseif mod(i,4) == 3
                    flag_dir = "top"
                    break;
                elseif mod(i,4) == 0
                    flag_dir = "right"
                    break;
                else
                    error("error2")
                end
            end 
        end
        return flag_end,flag_dir,r,x
    end    
    counter = 0;

    while flag_end == false && counter < 1e4
        if flag_dir == "top"
            flag_end,flag_dir,r,x = next_ind(2,r,x)
        elseif flag_dir == "left"
            flag_end,flag_dir,r,x = next_ind(1,r,x)
        elseif flag_dir == "bottom"
            flag_end,flag_dir,r,x = next_ind(4,r,x)
        elseif flag_dir == "right"
            flag_end,flag_dir,r,x = next_ind(3,r,x)
        else
            error("error")
        end
        counter = counter + 1
    end
    if counter == 1e4
        error("Counter reached 1e4")
    end
    return r,x
end


plot(scatter(x=r1,y=x1))

using FileIO
tmp = load("droop_I090_R_64_1_0_X_40_1_0_v2.jld")
tmp = load("droop_I100_R_111_1_0_X_71_1_0_v2.jld")
tmp = load("matching_I100_R_136_1_0_X_88_1_0_v2.jld")

Rv= tmp["Rverlauf"];
Xv = tmp["Xverlauf"];
XRv= tmp["XR"];
XRt= tmp["XR_t"];


r1,x1 = DetermineBoundary(XRv,Rv,Xv);
r2,x2 = plotxkrit(reverse(XRv),collect(reverse(Rv)),collect(reverse(Xv)));
plot([scatter(x=r1,y=x1),scatter(x=r2,y=x2)])

plot(surface(x=Rv,y=Xv,z=XRv))