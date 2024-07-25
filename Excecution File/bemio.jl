include("ReadWriteFILE/readWAMIT.jl")
include("functions/radiationIRF.jl")
include("functions/radiationIRFSS.jl")
include("functions/excitationIRF.jl")
include("functions/normalizeBEM.jl")
#include(ReadWriteFILE/writeBEMIOH5.jl")

## Test for the two case studies: RM3 and OSWEC. 
## option 0 is for RM3 case, option 1 is for OSWEC case
## Until the moment the only files not able to run correctly are the radiationIRFSS and writeBEMIOH5 functions


# User input 
option = 0

if option == 0
    # WAMIT
    hydro = readWAMIT("C:/Users/jelope/Desktop/Git/WEC-Sim/examples/RM3/hydroData/rm3.out")
    hydro = radiationIRF(hydro, 60)
    #hydro = radiationIRFSS(hydro)
    excitationIRF(hydro, 157)
    
    # Print field values for comparison
    print_hydro_data(hydro)

    #writeBEMIOH5(hydro)
    # Plot hydro data
    # plotBEMIO(hydro)
    
elseif option == 1
    # OSWEC
    hydro = readWAMIT("C:/Users/jelope/Desktop/Git/WEC-Sim/examples/OSWEC/hydroData/oswec.out")

    radiationIRF(hydro, 30)
    #radiationIRFSS(hydro)
    excitationIRF(hydro, 30)
    #writeBEMIOH5(hydro)

    print_hydro_data(hydro)

    # Plot hydro data
    # plotBEMIO(hydro)
else
    println("Invalid option selected. Please run the program again and select either 0 or 1.")
end
