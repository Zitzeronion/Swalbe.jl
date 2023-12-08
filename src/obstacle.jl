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