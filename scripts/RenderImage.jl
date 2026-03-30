using Images, FileIO

function data2grayscale(R, r, a; path = "data/")
    if R == 180 && r == 40 && a == 20
        id = 202531205
    elseif R == 120 && r == 80 && a == 20
        id = 20253112041
    end
    simset = (R, r, a, id)

    # hlist = FileIO.load("$(path)Rivulets/slip_non_qube_delta_05_height_R_$(simset[1])_r_$(simset[2])_ang_$(simset[3])_kbt_0.0_nm_3-2_runDate_$(simset[4]).jld2")
    hlist = FileIO.load(
        "$(path)Rivulets/slip_non_qube_height_R_$(simset[1])_r_$(simset[2])_ang_$(simset[3])_kbt_0.0_nm_3-2_runDate_$(simset[4]).jld2",
    )

    hmax = 0

    for i = 25000:25000:5000000
        hfield = hlist["h_$(i)"]
        if maximum(hfield) > hmax
            hmax = maximum(hfield)
        end
    end

    for i in enumerate(25000:25000:5000000)
        hf = reshape(hlist["h_$(i[2])"], 512, 512)
        hfnorm = hf/hmax
        grayH = Gray.(hfnorm)
        if i[1] < 10
            tlabel = "00$(i[1])"
        elseif (i[1] >= 10) && (i[1] < 100)
            tlabel = "0$(i[1])"
        else
            tlabel = "$(i[1])"
        end
        Images.save(
            "$(path)grayscale/R$(simset[1])_r$(simset[2])_a$(simset[3])/img_$(tlabel).png",
            grayH,
        )
    end
end
