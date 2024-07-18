using ProgressMeter, HDF5

function writeBEMIOH5(hydro::HydroData)
    # Writes the hydro data structure to a .h5 file.
    #
    # See ``WEC-Sim\tutorials\BEMIO`` for examples of usage.
    #
    # Parameters
    # ----------
    #     hydro : [1 x 1] struct
    #         Structure of hydro data that is written to ``hydro.file``

    @showprogress "Writing data in h5 format..." for i in 1:sum(hydro.Nb)
        N = 1 + (1 + 1 + 1) * hydro.Nb  # Rough division of tasks
        filename = hydro.file * ".h5"

        # Create h5 output file
        fileID = h5open(filename, "w")

        # Write string data
        writeH5Text(fileID, ["bem_data"], "code", hydro.code)
        for i in 1:hydro.Nb
            writeH5Text(fileID, ["body" * string(i), "properties"], "name", hydro.body[i])
        end
        if hydro.h == Inf
            writeH5Text(fileID, ["simulation_parameters"], "water_depth", "infinite") # A residual of the python based code
        end
        close(fileID)

        # Write body-independent wave parameters
        writeH5Parameter(filename, "/simulation_parameters/scaled", 0, "", "") # use >>body(#).bemioFlag = 0; in input
        writeH5Parameter(filename, "/simulation_parameters/g", hydro.g, "Gravitational acceleration", "m/s^2")
        writeH5Parameter(filename, "/simulation_parameters/rho", hydro.rho, "Water density", "kg/m^3")
        writeH5Parameter(filename, "/simulation_parameters/T", hydro.T, "Wave periods", "s")
        writeH5Parameter(filename, "/simulation_parameters/w", hydro.w, "Wave frequencies", "rad/s")
        if hydro.h != Inf     # A residual of the python based code
            writeH5Parameter(filename, "/simulation_parameters/water_depth", hydro.h, "Water depth", "m")
        end
        writeH5Parameter(filename, "/simulation_parameters/wave_dir", hydro.theta, "Wave direction", "deg")

        n = 0
        m_add = 0
        for i in 1:hydro.Nb
            m = hydro.dof[i]
        
            # Write body properties
            writeH5Parameter(filename, "/body" * string(i) * "/properties/body_number", i, "Number of rigid body from the BEM simulationCenter of gravity", "m")
            writeH5Parameter(filename, "/body" * string(i) * "/properties/cb", hydro.cb[:,i]', "Center of buoyancy", "m")
            writeH5Parameter(filename, "/body" * string(i) * "/properties/cg", hydro.cg[:,i]', "Center of gravity", "m")
            writeH5Parameter(filename, "/body" * string(i) * "/properties/disp_vol", hydro.Vo[i], "Displaced volume", "m^3")
            writeH5Parameter(filename, "/body" * string(i) * "/properties/dof", hydro.dof[i], "Degrees of freedom", "")
            writeH5Parameter(filename, "/body" * string(i) * "/properties/dof_start", m_add + 1, "Degrees of freedom", "")
            writeH5Parameter(filename, "/body" * string(i) * "/properties/dof_end", m_add + m, "Degrees of freedom", "")

            # Write GBM coefficients if available
            if haskey(hydro, :gbm)
                writeH5Parameter(filename, "/body" * string(i) * "/properties/mass", permutedims(hydro.gbm[(n+1):(n+m),:,1], (3, 2, 1)), "Generalized body modes mass", "kg")
                writeH5Parameter(filename, "/body" * string(i) * "/properties/damping", permutedims(hydro.gbm[(n+1):(n+m),:,2], (3, 2, 1)), "Generalized body modes damping", "N/m")
                writeH5Parameter(filename, "/body" * string(i) * "/properties/stiffness", permutedims(hydro.gbm[(n+1):(n+m),:,3], (3, 2, 1)), "Generalized body modes stiffness", "N-s/m")
            end

            # Write hydrostatic stiffness
            if haskey(hydro, :gbm)
                tmp = permutedims(hydro.gbm[(n+1):(n+m),:,4], (3, 2, 1))
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/linear_restoring_stiffness", tmp[1, m_add + 1:m_add + m, :], "Hydrostatic stiffness", "N/m")
                tmp = nothing
            else
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/linear_restoring_stiffness", hydro.Khs[:,:,i]', "Hydrostatic stiffness matrix", "")        
            end

            # Write added mass coefficients
            m_add += m
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/added_mass/inf_freq", permutedims(hydro.Ainf[(n+1):(n+m),:], (2, 1)), "Infinite frequency added mass", "kg")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/added_mass/all", permutedims(hydro.A[(n+1):(n+m),:,:], (3, 2, 1)), "Added mass", "kg-m^2 (rotation); kg (translation)")

            # Write excitation coefficients and IRF
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/im", permutedims(hydro.ex_im[(n+1):(n+m),:,:], (3, 2, 1)), "Imaginary component of excitation force", "")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/mag", permutedims(hydro.ex_ma[(n+1):(n+m),:,:], (3, 2, 1)), "Magnitude of excitation force", "")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/phase", permutedims(hydro.ex_ph[(n+1):(n+m),:,:], (3, 2, 1)), "Phase angle of excitation force", "rad")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/re", permutedims(hydro.ex_re[(n+1):(n+m),:,:], (3, 2, 1)), "Real component of excitation force", "")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/scattering/im", permutedims(hydro.sc_im[(n+1):(n+m),:,:], (3, 2, 1)), "Imaginary component of scattering force", "")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/scattering/mag", permutedims(hydro.sc_ma[(n+1):(n+m),:,:], (3, 2, 1)), "Magnitude of scattering force", "")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/scattering/phase", permutedims(hydro.sc_ph[(n+1):(n+m),:,:], (3, 2, 1)), "Phase angle of scattering force", "rad")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/scattering/re", permutedims(hydro.sc_re[(n+1):(n+m),:,:], (3, 2, 1)), "Real component of scattering force", "")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/froude-krylov/im", permutedims(hydro.fk_im[(n+1):(n+m),:,:], (3, 2, 1)), "Imaginary component of Froude-Krylov force", "")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/froude-krylov/mag", permutedims(hydro.fk_ma[(n+1):(n+m),:,:], (3, 2, 1)), "Magnitude of Froude-Krylov force", "")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/froude-krylov/phase", permutedims(hydro.fk_ph[(n+1):(n+m),:,:], (3, 2, 1)), "Phase angle of Froude-Krylov force", "rad")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/froude-krylov/re", permutedims(hydro.fk_re[(n+1):(n+m),:,:], (3, 2, 1)), "Real component of Froude-Krylov force", "")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/impulse_response_fun/f", permutedims(hydro.ex_K[(n+1):(n+m),:,:], (3, 2, 1)), "Impulse response function", "")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/impulse_response_fun/t", hydro.ex_t, "Time vector for the impulse response function", "s")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/impulse_response_fun/w", hydro.ex_w, "Interpolated frequencies used to compute the impulse response function", "s")

            # Write mean drift coefficients if available
            if haskey(hydro, :md_mc) # Only if mean drift variables (momentum conservation) have been calculated in BEM
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/mean_drift/momentum_conservation/val", permutedims(hydro.md_mc[(n+1):(n+m),:,:], (3, 2, 1)), "Value of mean drift force (momentum conservation)", "")
            end
            if haskey(hydro, :md_cs) # Only if mean drift variables (control surface approach) have been calculated in BEM
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/mean_drift/control_surface/val", permutedims(hydro.md_cs[(n+1):(n+m),:,:], (3, 2, 1)), "Value of mean drift force (control surface)", "")
            end
            if haskey(hydro, :md_pi) # Only if mean drift variables (pressure integration approach) have been calculated in BEM
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/mean_drift/pressure_integration/val", permutedims(hydro.md_pi[(n+1):(n+m),:,:], (3, 2, 1)), "Value of mean drift force (pressure integration)", "")
            end

            # Write radiation damping coefficients and IRF
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/all", permutedims(hydro.B[(n+1):(n+m),:,:], (3, 2, 1)), "Radiation damping", "")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/impulse_response_fun/K", permutedims(hydro.ra_K[(n+1):(n+m),:,:], (3, 2, 1)), "Impulse response function", "")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/impulse_response_fun/t", hydro.ra_t, "Time vector for the impulse response function", "s")
            writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/impulse_response_fun/w", hydro.ra_w, "Interpolated frequencies used to compute the impulse response function", "s")

            # Write radiation damping state space coefficients (optional)
            if haskey(hydro, :ss_A)
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/state_space/A/all", permutedims(hydro.ss_A[(n+1):(n+m),:,:,:], (4, 3, 2, 1)), "State Space A Coefficient", "")
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/state_space/B/all", permutedims(hydro.ss_B[(n+1):(n+m),:,:,:], (4, 3, 2, 1)), "State Space B Coefficient", "")
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/state_space/C/all", permutedims(hydro.ss_C[(n+1):(n+m),:,:,:], (4, 3, 2, 1)), "State Space C Coefficient", "")
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/state_space/D/all", permutedims(hydro.ss_D[(n+1):(n+m),:,:], (2, 1)), "State Space D Coefficient", "")
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/state_space/K/all", permutedims(hydro.ss_K[(n+1):(n+m),:,:,:], (4, 3, 2, 1)), "State Space K Coefficient", "")
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/state_space/it", permutedims(hydro.ss_O[(n+1):(n+m),:], (2, 1)), "Order of state space realization", "")
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/state_space/r2t", permutedims(hydro.ss_R2[(n+1):(n+m),:], (2, 1)), "State space curve fitting R**2 value", "")
                writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/state_space/conv", permutedims(hydro.ss_conv[(n+1):(n+m),:], (2, 1)), "State space conversion status", "")
            end

            for j in (n+1):(n+m)
                for k in 1:sum(hydro.dof)
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/added_mass/components/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.A[j,k,:], (3, 2, 1))]', "Added mass components as a function of frequency", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/components/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.B[j,k,:], (3, 2, 1))]', "Radiation damping components as a function of frequency", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/impulse_response_fun/components/K/" * string(j-m*i+m) * "_" * string(k), [hydro.ra_t', permutedims(hydro.ra_K[j,k,:], (3, 2, 1))]', "Components of the IRF", "")
                    if haskey(hydro, :ss_A) # Only if state space variables have been calculated
                        writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/state_space/A/components/" * string(j-m*i+m) * "_" * string(k), permutedims(hydro.ss_A[j,k,:,:], (4, 3, 2, 1)), "Components of the State Space A Coefficient", "")
                        writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/state_space/B/components/" * string(j-m*i+m) * "_" * string(k), permutedims(hydro.ss_B[j,k,:,:], (4, 3, 2, 1)), "Components of the State Space B Coefficient", "")
                        writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/state_space/C/components/" * string(j-m*i+m) * "_" * string(k), permutedims(hydro.ss_C[j,k,:,:], (4, 3, 2, 1)), "Components of the State Space C Coefficient", "")
                        writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/radiation_damping/state_space/D/components/" * string(j-m*i+m) * "_" * string(k), permutedims(hydro.ss_D[j,k], (2, 1)), "Components of the State Space D Coefficient", "")
                    end
                end
            end

            for j in (n+1):(n+m)
                for k in 1:hydro.Nh
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/components/im/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.ex_im[j,k,:], (3, 2, 1))]', "Imaginary component of excitation force as a function of frequency", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/components/re/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.ex_re[j,k,:], (3, 2, 1))]', "Real component of excitation force as a function of frequency", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/components/mag/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.ex_ma[j,k,:], (3, 2, 1))]', "Magnitude of excitation force as a function of frequency", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/components/phase/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.ex_ph[j,k,:], (3, 2, 1))]', "Phase of excitation force as a function of frequency", "")
                    if haskey(hydro, :md_mc) # Only if mean drift variables (momentum conservation) have been calculated in BEM
                        writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/mean_drift/momentum_conservation/components/val/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.md_mc[j,k,:], (3, 2, 1))]', "Magnitude of mean drift force (momentum conservation) as a function of frequency", "")
                    end
                    if haskey(hydro, :md_cs) # Only if mean drift variables (control surface approach) have been calculated in BEM
                        writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/mean_drift/control_surface/components/val/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.md_cs[j,k,:], (3, 2, 1))]', "Magnitude of mean drift force (control surface) as a function of frequency", "")
                    end
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/impulse_response_fun/components/f/" * string(j-m*i+m) * "_" * string(k), [hydro.ex_t', permutedims(hydro.ex_K[j,k,:], (3, 2, 1))]', "Components of the IRF:f", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/scattering/components/im/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.sc_im[j,k,:], (3, 2, 1))]', "Imaginary component of scattering force as a function of frequency", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/scattering/components/mag/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.sc_ma[j,k,:], (3, 2, 1))]', "Magnitude of scattering force as a function of frequency", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/scattering/components/phase/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.sc_ph[j,k,:], (3, 2, 1))]', "Phase of scattering force as a function of frequency", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/scattering/components/re/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.sc_re[j,k,:], (3, 2, 1))]', "Real component of scattering force as a function of frequency", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/froude-krylov/components/im/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.fk_im[j,k,:], (3, 2, 1))]', "Imaginary component of Froude-Krylov force as a function of frequency", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/froude-krylov/components/mag/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.fk_ma[j,k,:], (3, 2, 1))]', "Magnitude of Froude-Krylov force as a function of frequency", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/froude-krylov/components/phase/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.fk_ph[j,k,:], (3, 2, 1))]', "Phase of Froude-Krylov force as a function of frequency", "")
                    writeH5Parameter(filename, "/body" * string(i) * "/hydro_coeffs/excitation/froude-krylov/components/re/" * string(j-m*i+m) * "_" * string(k), [hydro.T', permutedims(hydro.fk_re[j,k,:], (3, 2, 1))]', "Real component of Froude-Krylov force as a function of frequency", "")
                end
            end

            n += m
        end
    end
end

function writeH5Parameter(filename, location, values, description, units)
    h5write(filename, location, values)
    h5writeattr(filename, location, "description", description)
    h5writeattr(filename, location, "units", units)
end

function writeH5Text(fileID, group, dataset, data)
    if length(group) == 1
        group_id = H5G.create(fileID, group[1], "H5P_DEFAULT", "H5P_DEFAULT", "H5P_DEFAULT")
        space_id = H5S.create("H5S_SCALAR")
        stype = H5T.copy("H5T_C_S1")
        H5T.set_size(stype, length(data))
        dataset_id = H5D.create(group_id, dataset, stype, space_id, "H5P_DEFAULT")
        H5D.write(dataset_id, stype, "H5S_ALL", "H5S_ALL", "H5P_DEFAULT", data)
        H5D.close(dataset_id)
        H5S.close(space_id)
        H5G.close(group_id)
    elseif length(group) == 2
        group_1_id = H5G.create(fileID, group[1], "H5P_DEFAULT", "H5P_DEFAULT", "H5P_DEFAULT")
        group_2_id = H5G.create(group_1_id, group[2], "H5P_DEFAULT", "H5P_DEFAULT", "H5P_DEFAULT")
        space_id = H5S.create("H5S_SCALAR")
        stype = H5T.copy("H5T_C_S1")
        H5T.set_size(stype, length(data))
        dataset_id = H5D.create(group_2_id, dataset, stype, space_id, "H5P_DEFAULT")
        H5D.write(dataset_id, stype, "H5S_ALL", "H5S_ALL", "H5P_DEFAULT", data)
        H5D.close(dataset_id)
        H5S.close(space_id)
        H5G.close(group_2_id)
        H5G.close(group_1_id)
    end
end
