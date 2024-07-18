function normalizeBEM(data::HydroData)
    # Sort, if necessary
    if !issorted(data.w)
        I = sortperm(data.w)
        data.w = data.w[I]
        data.T = data.T[I]
        data.A = data.A[:, :, I]
        data.B = data.B[:, :, I]
        data.ex_ma = data.ex_ma[:, :, I]
        data.ex_ph = data.ex_ph[:, :, I]
        data.ex_re = data.ex_re[:, :, I]
        data.ex_im = data.ex_im[:, :, I]
        data.sc_ma = data.sc_ma[:, :, I]
        data.sc_ph = data.sc_ph[:, :, I]
        data.sc_re = data.sc_re[:, :, I]
        data.sc_im = data.sc_im[:, :, I]
        data.fk_ma = data.fk_ma[:, :, I]
        data.fk_ph = data.fk_ph[:, :, I]
        data.fk_re = data.fk_re[:, :, I]
        data.fk_im = data.fk_im[:, :, I]
    end

    if data.code != "WAMIT"  # Normalize
        g = data.g
        rho = data.rho
        data.Khs = data.Khs / (g * rho)
        data.A = data.A / rho
        data.Ainf = data.A[:, :, end]  # Assuming this is correct based on your structure
        for i in 1:length(data.w)
            data.B[:, :, i] = data.B[:, :, i] / (rho * data.w[i])
        end
        data.ex_ma = data.ex_ma / (g * rho)
        data.ex_re = data.ex_re / (g * rho)
        data.ex_im = data.ex_im / (g * rho)
        data.sc_ma = data.sc_ma / (g * rho)
        data.sc_re = data.sc_re / (g * rho)
        data.sc_im = data.sc_im / (g * rho)
        data.fk_ma = data.fk_ma / (g * rho)
        data.fk_re = data.fk_re / (g * rho)
        data.fk_im = data.fk_im / (g * rho)
    end

    return data
end
