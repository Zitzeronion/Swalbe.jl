
# """
# function obslist(obs)
# Returns a list containing the relevgant obstacle nodes constructed from interior, obsleft, obsright, obsup, obsdown, corneroutlu, corneroutld, corneroutru, corneroutrd, cornerru, cornerrd, cornerlu, cornerld

# TODO: Find out how double corners e.g. leftup and leftdown should be handled 
#     Add a field for outer corners. 
# """
# function obslist(obs; T=Float64, verbose=true)
#         if verbose 
#             println("WARNING: Always make your obstacles two nodes thick or the algorithm will crash.")
#         end
#         Lx = size(obs)[1]
#         Ly = size(obs)[2]
#         interior::Array{T,2} = zeros(Lx,Ly)
#         obsleft::Array{T,2} = zeros(Lx, Ly)
#         obsright::Array{T,2} = zeros(Lx, Ly)
#         obsup::Array{T,2} = zeros(Lx, Ly)
#         obsdown::Array{T,2} = zeros(Lx, Ly)
#         corneroutru::Array{T,2} = zeros(Lx, Ly)
#         corneroutrd::Array{T,2} = zeros(Lx, Ly)
#         corneroutlu::Array{T,2} = zeros(Lx, Ly)
#         corneroutld::Array{T,2} = zeros(Lx, Ly)
#         cornerru::Array{T,2} = zeros(Lx, Ly)
#         cornerrd::Array{T,2} = zeros(Lx, Ly)
#         cornerlu::Array{T,2} = zeros(Lx, Ly)
#         cornerld::Array{T,2} = zeros(Lx, Ly)

#         # The names where given as follows

#         #
#         #  For a empty box
#         #        xxxxxxxx  obsup
#         #        x      x
#         # obsleftx      x obsright
#         #        xxxxxxxx obsdown
#         #  
#         # Which makes the names weird for a filled shape
#         #           obsdown
#         #          xxxxxxxx  
#         # obsright xxxxxxxx obsleft
#         #          xxxxxxxx 
#         #          xxxxxxxx 
#         #           obsup
#         #Corners are named like this
#         # corneroutlu      corneroutru
#         #          xxxxxxxx  
#         #          xxxxxxxx
#         #          xxxxxxxx 
#         #          xxxxxxxx 
#         # corneroutld      corneroutrd
#         #
#         # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#         # x cornerlu              cornerru  x 
#         # x                                 X
#         # x cornerld               cornerrd x 
#         # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#         for i in 1:Lx, j in 1:Ly
#             interior[i,j] = (obs[i,j]==1 && obs[i,mod1(j-1,Ly)]==1  && obs[i,mod1(j+1,Ly)]==1 && obs[mod1(i-1,Lx),j]==1 && obs[mod1(i-1,Lx),mod1(j-1,Ly)]==1 && obs[mod1(i-1,Lx),mod1(j+1,Ly)]==1 && obs[mod1(i+1,Lx),j]==1 && obs[mod1(i+1,Lx),mod1(j-1,Ly)]==1 && obs[mod1(i+1,Lx),mod1(j+1,Ly)]==1) ? 1 : 0 
#             #-------------------------------------Linesegment-----------------------------------------corner----------------------------------left up corner-------------------------------------------------------------------------------left down corner--------------------------------------------------
#             obsleft[i,j] = (obs[i,j]==1 && ((obs[i, mod1(j+1,Ly)]==0) || (obs[i, mod1(j+1,Ly)]==1 && (  (obs[mod1(i+1,Lx),mod1(j+1,Ly)]==0 && obs[mod1(i+1,Lx),j]==1)   || (obs[mod1(i-1,Lx),mod1(j+1,Ly)]==0 && obs[mod1(i-1,Lx),j]==1) )))) ? 1 : 0
#             #-------------------------------------Linesegment-----------------------------------------corner----------------------------------right up corner-------------------------------------------------------------------------------right down corner--------------------------------------------------
#             obsright[i,j] = (obs[i,j]==1 && ((obs[i, mod1(j-1,Ly)]==0) || (obs[i, mod1(j-1,Ly)]==1 && (  (obs[mod1(i+1,Lx),mod1(j-1,Ly)]==0 && obs[mod1(i+1,Lx),j]==1)   || (obs[mod1(i-1,Lx),mod1(j-1,Ly)]==0 && obs[mod1(i-1,Lx),j]==1) )))) ? 1 : 0
#             #-------------------------------------Linesegment-----------------------------------------corner----------------------------------right up corner-------------------------------------------------------------------------------left up corner--------------------------------------------------
#             obsup[i,j] = (obs[i,j]==1 && ((obs[mod1(i+1,Lx), j]==0) || (obs[mod1(i+1,Lx), j]==1 && (  (obs[mod1(i+1,Lx),mod1(j-1,Ly)]==0 && obs[i,mod1(j-1,Ly)]==1)   || (obs[mod1(i+1,Lx),mod1(j+1,Ly)]==0 && obs[i,mod1(j+1,Ly)]==1)  )))) ? 1 : 0
#             #-------------------------------------Linesegment-----------------------------------------corner----------------------------------right down corner-----------------------------------------------------------------------------left down corner--------------------------------------------------
#             obsdown[i,j] = (obs[i,j]==1 && ((obs[mod1(i-1,Lx), j]==0) || (obs[mod1(i-1,Lx), j]==1 && (  (obs[mod1(i-1,Lx),mod1(j-1,Ly)]==0 && obs[i,mod1(j-1,Ly)]==1)   || (obs[mod1(i-1,Lx),mod1(j+1,Ly)]==0 && obs[i,mod1(j+1,Ly)]==1)  )))) ? 1 : 0
#             #outer corners
#             corneroutlu[i,j] = ( obs[i,j]==1 && obs[mod1(i-1,Lx),j]==0 &&  obs[mod1(i-1,Lx),mod1(j-1,Ly)]==0 && obs[i,mod1(j-1,Ly)]==0) ? 1 : 0
#             corneroutld[i,j] = ( obs[i,j]==1 && obs[mod1(i+1,Lx),j]==0 &&  obs[mod1(i+1,Lx),mod1(j-1,Ly)]==0 && obs[i,mod1(j-1,Ly)]==0) ? 1 : 0
#             corneroutru[i,j] = ( obs[i,j]==1 && obs[mod1(i-1,Lx),j]==0 &&  obs[mod1(i-1,Lx),mod1(j+1,Ly)]==0 && obs[i,mod1(j+1,Ly)]==0) ? 1 : 0
#             corneroutrd[i,j] = ( obs[i,j]==1 && obs[mod1(i+1,Lx),j]==0 &&  obs[mod1(i+1,Lx),mod1(j+1,Ly)]==0 && obs[i,mod1(j+1,Ly)]==0) ? 1 : 0
#         end
#         #inner corners
#         cornerlu = obsleft .* obsup .- corneroutrd
#         cornerld = obsleft .* obsdown .- corneroutru
#         cornerru = obsright .* obsup .- corneroutld
#         cornerrd = obsright .* obsdown .- corneroutlu


#         #preparing return
#         border1 = zeros(Lx, Ly)
#         border2 = zeros(Lx, Ly)
#         border3 = zeros(Lx, Ly)
#         border4 = zeros(Lx, Ly)
#         border5 = zeros(Lx, Ly)
#         border6 = zeros(Lx, Ly)
#         border7 = zeros(Lx, Ly)
#         border8 = zeros(Lx, Ly)     

#         border1 .= obsright .- corneroutld .- corneroutlu .+ interior
#         border2 .= obsup .- corneroutld .- corneroutrd .+ interior
#         border3 .= obsleft .- corneroutrd .- corneroutru .+ interior
#         border4 .= obsdown .- corneroutru .- corneroutlu .+ interior
#         border5 .= obsright .+ obsup .- cornerru .- corneroutld .- corneroutrd .- corneroutlu .+ interior
#         border6 .= obsleft .+ obsup .- cornerlu .- corneroutld .- corneroutrd .- corneroutru .+ interior
#         border7 .= obsleft .+ obsdown .- cornerld .- corneroutlu .- corneroutrd .- corneroutru .+ interior
#         border8 .= obsright .+ obsdown .- cornerrd .- corneroutld .- corneroutlu .- corneroutru .+ interior

#         #retrun  interior, obsleft, obsright, obsup, obsdown, corneroutlu, corneroutld, corneroutru, corneroutrd, cornerru, cornerrd, cornerlu, cornerld, [border1, border2, border3, border4, border5, border6, border7, border8] 
#         return interior, [border1, border2, border3, border4, border5, border6, border7, border8] 
# end



# """
#     function obslist!(sys::SysConstWithBound; T=Float64, verbose=false)
# Updates interior and border according to obs. See obslist().
# """
# function obslist!(sys::SysConstWithBound; T=Float64, verbose=false)
#     list = obslist(sys.obs, T=T, verbose=verbose)
#     sys.interior .= list[1]
#     sys.border[1] .= list[2][1]
#     sys.border[2] .= list[2][2]
#     sys.border[3] .= list[2][3]
#     sys.border[4] .= list[2][4]
#     sys.border[5] .= list[2][5]
#     sys.border[6] .= list[2][6]
#     sys.border[7] .= list[2][7]
#     sys.border[8] .= list[2][8]


#       # Straight elements j+1, i+1, i-1, j-1
#       oip = circshift(sys.interior, (1,0))
#       ojp = circshift(sys.interior, (0,1))
#       oim = circshift(sys.interior, (-1,0))
#       ojm = circshift(sys.interior, (0,-1))
#       # Diagonal elements  
#       oipjp = circshift(sys.interior, (1,1))
#       oimjp = circshift(sys.interior, (-1,1))
#       oimjm = circshift(sys.interior, (-1,-1))
#       oipjm = circshift(sys.interior, (1,-1))
#       #One could think about putting weights here that are less for the diagonal ones than for the straight ones
#     sys.solidneighbours .= 1 ./ max.(1 .- 1/8 .* (oip .+ ojp .+ oim .+ ojm .+ oipjp .+ oimjp .+ oimjm .+ oipjm), 1/8)
#     return nothing
# end




"""
function obslist1D(obs)
Returns a list containing the relevgant obstacle nodes constructed from interior, obsleft, obsright, obsup, obsdown, corneroutlu, corneroutld, corneroutru, corneroutrd, cornerru, cornerrd, cornerlu, cornerld

"""
function obslist1D(obs; T=Float64, verbose=true)
        if verbose 
            println("WARNING: Always make your obstacles four nodes thick or the algorithm will crash.")
        end
        L = size(obs)[1]
        interior::Array{T} = zeros(L)
        obsleft::Array{T} =  zeros(L)
        obsright::Array{T} = zeros(L)

         # The names where given as follows
        #
        #  For a empty box
        # obsleft xx      xx obsright
        #  
        # Which makes the names weird for a filled shape
        # obsright xxxxxxxx obsleft

        for i in 1:L 
            interior[i] = (obs[i]==1 && obs[mod1(i+1,L)]==1 && obs[mod1(i-1,L)]==1) ? 1 : 0
            obsleft[i] = (obs[i]==1 && obs[mod1(i-1,L)]== 1) ? 1 : 0
            obsright[i] =  (obs[i]==1 && obs[mod1(i+1,L)]== 1) ? 1 : 0        
        end


      
        return interior, [obsright, obsleft] 
end



"""
    function obslist1D!(sys::SysConstWithBound; T=Float64, verbose=false)
Updates interior and border according to obs. See obslist().
"""
function obslist!(sys::SysConstWithBound_1D; T=Float64, verbose=false)
    list = obslist1D(sys.obs, T=T, verbose=verbose)
    sys.interior .= list[1]
    sys.border[1] .= list[2][1]
    sys.border[2] .= list[2][2]


      # Straight elements j+1, i+1, i-1, j-1
      oip = circshift(sys.interior, 1)
      oim = circshift(sys.interior, -1)
      #One could think about putting weights here that are less for the diagonal ones than for the straight ones
    return nothing
end



# """
#     function obstacleheight!(state::StateWithBound, sys:: SysConstWithBound)
# Gives the height a value inside the obstacle that is calculated with respect to the contact angle of the substrate of the obstacle surface. 
# When calculating the filmpressure this not physically existent film height inside obstacle nodes will via Laplacian of height incluence the 
# surface tension and thgerefore create the miniscus obsverable near a wall like obstacle. 

# # Mathematics

# TODO: 

# # ATENTION
# Induces mass to the system
# """
# function obstacleheight!(state::StateWithBound, sys:: SysConstWithBound)
#       hip, hjp, him, hjm, hipjp, himjp, himjm, hipjm = Swalbe.viewneighbors(state.dgrad)
#       oip, ojp, oim, ojm, oipjp, oimjp, oimjm, oipjm = Swalbe.viewneighbors(state.obsgrad)
#       # Straight elements j+1, i+1, i-1, j-1
#       hip = circshift(state.height, (1,0))
#       hjp = circshift(state.height, (0,1))
#       him = circshift(state.height, (-1,0))
#       hjm = circshift(state.height, (0,-1))
#       # Diagonal elements  
#       hipjp = circshift(state.height, (1,1))
#       himjp = circshift(state.height, (-1,1))
#       himjm = circshift(state.height, (-1,-1))
#       hipjm = circshift(state.height, (1,-1))
#         # Straight elements j+1, i+1, i-1, j-1
#       oip = circshift(sys.interior, (1,0))
#       ojp = circshift(sys.interior, (0,1))
#       oim = circshift(sys.interior, (-1,0))
#       ojm = circshift(sys.interior, (0,-1))
#       # Diagonal elements  
#       oipjp = circshift(sys.interior, (1,1))
#       oimjp = circshift(sys.interior, (-1,1))
#       oimjm = circshift(sys.interior, (-1,-1))
#       oipjm = circshift(sys.interior, (1,-1))

#       #Again other weights than 1/8 are possible see obslist!()
#       state.height .= ((1  .- sys.interior) .* state.height 
#                         .+ sys.interior .* ( sys.solidneighbours .* 1/8 .* ((1 .- oip) .* hip .+ (1 .- ojp) .* hjp .+ (1 .- oim) .* him .+ (1 .- ojm) .* hjm
#                                                                      .+ (1 .-oipjp) .* hipjp .+ (1 .- oimjp) .* himjp .+ (1 .- oimjm) .* himjm .+ (1 .- oipjm) .* hipjm )
#                         .- tan(sys.θ * pi)
#                         )
#                 )
                                                                     
#     return nothing
# end



# function obstaclepressure!(state::StateWithBound, sys:: SysConstWithBound)
#     hip, hjp, him, hjm, hipjp, himjp, himjm, hipjm = Swalbe.viewneighbors(state.dgrad)
#     oip, ojp, oim, ojm, oipjp, oimjp, oimjm, oipjm = Swalbe.viewneighbors(state.obsgrad)
#     # Straight elements j+1, i+1, i-1, j-1
#     hip = circshift(state.pressure, (1,0))
#     hjp = circshift(state.pressure, (0,1))
#     him = circshift(state.pressure, (-1,0))
#     hjm = circshift(state.pressure, (0,-1))
#     # Diagonal elements  
#     hipjp = circshift(state.pressure, (1,1))
#     himjp = circshift(state.pressure, (-1,1))
#     himjm = circshift(state.pressure, (-1,-1))
#     hipjm = circshift(state.pressure, (1,-1))
#       # Straight elements j+1, i+1, i-1, j-1
#     oip = circshift(sys.interior, (1,0))
#     ojp = circshift(sys.interior, (0,1))
#     oim = circshift(sys.interior, (-1,0))
#     ojm = circshift(sys.interior, (0,-1))
#     # Diagonal elements  
#     oipjp = circshift(sys.interior, (1,1))
#     oimjp = circshift(sys.interior, (-1,1))
#     oimjm = circshift(sys.interior, (-1,-1))
#     oipjm = circshift(sys.interior, (1,-1))

#     #Again other weights than 1/8 are possible see obslist!()
#     state.pressure .= ((1  .- sys.interior) .* state.pressure 
#                       .+ sys.interior .* ( sys.solidneighbours .* 1/8 .* ((1 .- oip) .* hip .+ (1 .- ojp) .* hjp .+ (1 .- oim) .* him .+ (1 .- ojm) .* hjm
#                                                                    .+ (1 .-oipjp) .* hipjp .+ (1 .- oimjp) .* himjp .+ (1 .- oimjm) .* himjm .+ (1 .- oipjm) .* hipjm )
#                       #.+ tan(sys.θ * pi) Maybe we could also add somethign. But adding stuff is what kills mass conservation so better don't
#                       )
#               )
                                                                   
#   return nothing
# end


# """
#     function tidyupobs!(state::StateWithBound, sys::SysConstWithBound)
# just sets constant height zero inside obstacle to have nicer looking data and not accidently analyse solid nodes data
# """
# function tidyupobs!(state::StateWithBound, sys::SysConstWithBound)
#     state.height .= (1 .- sys.obs) .* state.height .+ sys.obs .*2
# end