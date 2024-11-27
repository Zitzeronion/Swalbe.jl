"""
    view_four()

Splits a chuck of memory in four equivalent chucks
"""
function view_four(f)
    f1 = view(f, :, 1)
    f2 = view(f, :, 2)
    f3 = view(f, :, 3)
    f4 = view(f, :, 4)
    
    return f1, f2, f3, f4
end


"""
    power_broad(arg, n)

Computes `arg` to the power `n`.

Actually this is useful because the `^` operator is much slower.
Same thing I learned about the `pow` function in **C**, * yes it does what you want, but it is slow as fuck *.

# Examples
```jldoctest
julia> using Swalbe, Test

julia> Swalbe.power_broad(3, 3)
27

julia> Swalbe.power_broad.([2.0 5.0 6.0], 2) # Use the broadcasting operator `.`
1×3 Matrix{Float64}:
 4.0  25.0  36.0

```

See also: [`filmpressure`](@ref)
"""
function power_broad(arg::Float64, n::Int)
    temp = arg
    for i = 1:(n-1)
        temp *= arg
    end
    return temp
end

function power_broad(arg::Float32, n::Int)
    temp = arg
    for i = 1:(n-1)
        temp *= arg
    end
    return temp
end

function power_broad(arg::Int, n::Int)
    temp = arg
    for i = 1:(n-1)
        temp *= arg
    end
    return temp
end



"""
    power_2(arg::Float64)

Computes the square of `arg`, i.e., `arg^2`.

# Arguments
- `arg :: Float64`: The number to be squared.

# Returns
- The square of `arg`.

# Example
```jldoctest
julia> power_2(4)
16
```
"""
function power_2(arg::Float64)
    return arg * arg
end


"""
    power_3(arg::Float64)

Computes the cube of `arg`, i.e., `arg^3`.

# Arguments
- `arg :: Float64`: The number to be cubed.

# Returns
- The cube of `arg`.

# Example
```jldoctest
julia> power_3(2)
8
```
"""
function power_3(arg::Float64)
    return arg * arg * arg
end

"""
    power_4(arg::Float64)

Computes the fourth power of `arg`, i.e., `arg^4`. 

This implementation is approximately 0.25% faster than computing `power_2(power_2(arg))`, though the reason for this speedup is not immediately clear.

# Arguments
- `arg :: Float64`: The number to be raised to the fourth power.

# Returns
- The fourth power of `arg`.

# Example
```jldoctest
julia> power_4(2)
16
```
"""
function power_4(arg::Float64)
    return arg * arg * arg * arg
end


"""
    power_5(arg::Float64)

Computes the fifth power of `arg`, i.e., `arg^5`.

# Arguments
- `arg :: Float64`: The number to be raised to the fifth power.

# Returns
- The fifth power of `arg`.

# Example
```jldoctest
julia> power_5(2)
32
```
"""
function power_5(arg::Float64)
    return arg * arg * arg * arg * arg
end


"""
    power_6(arg::Float64)

Computes the sixth power of `arg` by using a combination of `power_2` and `power_3`, i.e., `arg^6 = power_2(power_3(arg))`.

# Arguments
- `arg :: Float64`: The number to be raised to the sixth power.

# Returns
- The sixth power of `arg`.

# Example
```jldoctest
julia> power_6(2)
64
```
"""
function power_6(arg::Float64)
    return power_2(power_3(arg))
end


"""
    power_8(arg::Float64)

Computes the eighth power of `arg`, i.e., `arg^8`. This is done by raising `arg` to the power of 2 four times, i.e., `arg^8 = power_2(power_2(power_2(arg)))`.

# Arguments
- `arg :: Float64`: The number to be raised to the eighth power.

# Returns
- The eighth power of `arg`.

# Example
```jldoctest
julia> power_8(2)
256
```
"""
function power_8(arg::Float64)
    return power_2(power_2(power_2(arg)))
end

"""
    power_12(arg::Float64)

Computes the twelfth power of `arg`, i.e., `arg^12`. This is done by combining `power_3` and `power_2` functions, i.e., `arg^12 = power_3(power_2(power_2(arg)))`.

# Arguments
- `arg :: Float64`: The number to be raised to the twelfth power.

# Returns
- The twelfth power of `arg`.

# Example
```jldoctest
julia> power_12(2)
4096
```
"""
function power_12(arg::Float64)
    return power_3(power_2(power_2(arg)))
end


"""
    power_9(arg::Float64)

Computes the ninth power of `arg`, i.e., `arg^9`. This is done by combining `power_3` function twice, i.e., `arg^9 = power_3(power_3(arg))`.

# Arguments
- `arg :: Float64`: The number to be raised to the ninth power.

# Returns
- The ninth power of `arg`.

# Example
```jldoctest
julia> power_9(2)
512
```
"""
function power_9(arg::Float64)
    return power_3(power_3(arg))
end

"""
    power_16(arg::Float64)

Computes the sixteenth power of `arg`, i.e., `arg^16`. This is done by applying the `power_4` function twice, i.e., `arg^16 = power_4(power_4(arg))`.

# Arguments
- `arg :: Float64`: The number to be raised to the sixteenth power.

# Returns
- The sixteenth power of `arg`.

# Example
```jldoctest
julia> power_16(2)
65536
```
"""
function power_16(arg::Float64)
    return power_4(power_4(arg))
end


"""
    fast_32(arg::Float64)

Computes `arg^3 - arg^2` more efficiently using the `power_3` and `power_2` functions. This method avoids computing `arg^2` multiple times.

# Arguments
- `arg :: Float64`: The input value.

# Returns
- The result of `arg^3 - arg^2`.

# Example
```jldoctest
julia> fast_32(2)
4
```
"""
function fast_32(arg::Float64)
    return power_3(arg) - power_2(arg)
end



"""
    fast_93(arg::Float64)

Computes `arg^3 - arg^3(arg)` more efficiently using the `power_3` function twice. This method avoids redundant power calculations.

# Arguments
- `arg :: Float64`: The input value.

# Returns
- The result of `arg^3 - arg^3(arg)`.

# Example
```jldoctest
julia> fast_93(2)
56
```
"""
function fast_93(arg::Float64)
    temp = power_3(arg)
    return power_3(temp) - temp
end


"""
    to_4_digits(i)

Converts an integer `i` into a 4-digit string, padding with leading zeros if necessary.

# Arguments
- `i :: Int`: The integer value to be formatted.

# Returns
- A string representation of `i` with at least 4 digits, padded with leading zeros if necessary.

# Example
```jldoctest
julia> to_4_digits(5)
"0005"

julia> to_4_digits(50)
"0050"

julia> to_4_digits(123)
"0123"
```
"""
function to_4_digits(i)
    if i < 10
        return "000$(i)"
    elseif i < 100
        return "00$(i)"
    elseif i < 1000
        return "0$(i)"
    else
        return "$(i)"
    end
end

function fast_63(arg::Float64, fac::Float64)
    temp = power_3(arg)
    # This was introducing anisotropy for reasons I don't understand
    # return fac*power_2(temp)-temp
    return power_2(temp)-temp
end

function wetting_potential_93(arg::Float64)
    temp = power_2(arg)
    return 8*temp - 2*power_4(temp)
end

function wetting_potential_63(arg::Float64)
    return 5*power_2(arg) - 2*power_5(arg)
end

function wetting_potential_32(arg::Float64)
    return 2*arg - arg*arg
    # return 0
end


"""
	function modified_modulo(x, y)

Calculates x%y for 1 base arrays, i.e. does not return 0 but y
"""
function modified_modulo(x, y)
	result = mod(x, y)
	return result == 0 ? y : result
end
