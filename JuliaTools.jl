"""
    Collection of functions for analysis of BRW (HDF5) files generated with the BrainWave program from the company 3Brain.
    Laboratory 19 of the CINVESTAV in charge of Dr. Rafael Gutierrez Aguilar.
    Work developed by Isabel Romero-Maldonado (2020 - )
    isabelrm.biofisica@gmail.com
    https://github.com/LBitn
"""

module JuliaTools

using HDF5
using JLD
using DelimitedFiles
using DataFrames
using Statistics
using Distributions
using BinningAnalysis
using Plots
using HistogramThresholding
using ImageContrastAdjustment
using Suppressor
using StatsBase
using SignalAnalysis
using DSP

export div_ab
export size_segments
export Get_chunks
export SavePaths
export VariablesBRW
export searchdir
export Digital2Analogue
export donoho
export Sparsity
export Density
export neighborgs
export MeanΔxCI
export Get_Groups
export FillingHolesCrux
export FigureGroups
export NuevePlots
export Thresholding
export Zplot
export ZW
export Z0
export ΔV
export SaturacionNegativaTemporal
export SegmentsComplet
export ms2frames
export FiltroMUAremez
export SaturacionPositiva
export TimesPotentialDischarge


# -------------------------------------------------------------------------------------------------------- #
"""
    searchdir( path::String, key::String ) -> Vector{String}
    Find inside the given path, the files with the key word on their name and returns a vector of strings with the full name of those files (complete path)
"""
@inline searchdir( path::String, key::String ) = filter( x -> endswith( x, key ), readdir( path; join = true ) );
# -------------------------------------------------------------------------------------------------------- #
# using HDF5
# using JLD
# using DelimitedFiles
"""
VariablesBRW( FILEBRW::String ) -> Variables::Dict{ String, Any }, FILEPATHS::String
    Reads the version of BRW and if it is the latest ( 102 ) version, extract a dictionary with the variables of interest
    It also creates an Info folder where it saves a .jld file with that dictionary for future uses.

    Notes: falta considerar las versiones viejas de brw, esto funciona para 2020 ->
"""
function VariablesBRW( FILEBRW::String )

    BRW = h5open( FILEBRW, "r" );
    Description = read( open_attribute( BRW, "Description" ) );
    BRWversion = read_attribute( open_group( BRW, "3BData" ), "Version" );
    BRWsize = ( ( stat( FILEBRW ).size ) / 1000000 ) / 1024;
    PATHBRWs = dirname( FILEBRW );
    if BRWversion == 102

        aux = split( Description, "," )[ end - 1 : end ];
        day = split( aux[ 1 ] )[ 2 ];
        month = split( aux[ 1 ] )[ 1 ][ 1 : 3 ];
        year = split( aux[ end ] )[ 1 ];
        brw = replace( splitpath( FILEBRW )[ end ], ".brw" => "" );
        PATHMain = joinpath( dirname( PATHBRWs ), string( brw, "-", day, "-", month, "-", year ) );
        mkpath( PATHMain );
        RecVars = read( open_group(  BRW, "3BRecInfo" ), "3BRecVars" );

        for i in keys( RecVars )
            RecVars[ i ] = RecVars[ i ][ 1 ];
        end

        MeaChip = read( open_group(  BRW, "3BRecInfo" ), "3BMeaChip" );
        delete!( MeaChip, "ROIs" ); delete!( MeaChip, "SysChs" );

        for i in keys( MeaChip )
            MeaChip[ i ] = MeaChip[ i ][ 1 ];
        end

        for i in keys( MeaChip )
            MeaChip[ i ] = convert( Int, MeaChip[ i ] );
        end

        MeaChip[ "Layout" ] = reshape( collect( 1:MeaChip[ "NRows" ] * MeaChip[ "NCols" ] ), MeaChip[ "NRows" ], MeaChip[ "NCols" ] );
        Noise = read( open_group( open_group(  BRW, "3BData" ), "3BInfo"), "3BNoise" );
        RawPath = Dict( "Raw" => "3BData/Raw" ); filename = h5open( FILEBRW, "r" ).filename;
        Description = Dict( "Description" => Description, "Version" => BRWversion, "BRW" => filename );
        RecVars[ "BitDepth" ] = convert( Int, RecVars[ "BitDepth" ] );
        SignalInversion = RecVars[ "SignalInversion" ];
        MaxVolt = RecVars[ "MaxVolt" ]; MinVolt = RecVars[ "MinVolt" ]; BitDepth = RecVars[ "BitDepth" ];
        ADCCountsToMV = SignalInversion * ( ( MaxVolt - MinVolt )/ 2^BitDepth );
        MVOffset = SignalInversion * MinVolt;
        Extras = Dict( "MVOffset" => MVOffset, "ADCCountsToMV" => ADCCountsToMV, "BRWsizeGB" => BRWsize );
        Variables = merge( Description, RecVars, MeaChip, Noise, RawPath, Extras );

        for i in keys( Variables )
            try
                Variables[ i ] = convert( Int, Variables[ i ] );
            catch
                Variables[ i ] = Variables[ i ];
            end
        end

        PATHInfo = joinpath( PATHMain, "Info" ); mkpath( PATHInfo );
        FILEVARS = joinpath( PATHInfo, "variablesBRW.jld");
        save( FILEVARS, "Variables", Variables );

        close( BRW )
        cd( PATHMain )

        println( "You are now working on the new main path: ", PATHMain );
        println( "with the file: ", basename( BRW.filename ) );
        println( Description[ "Description" ] );
        println( "HDF5 file size: $BRWsize GB" );

        TEXTVars = copy( Variables );
        delete!( TEXTVars, "StdMean" ); delete!( TEXTVars, "ValidChs" ); delete!( TEXTVars, "Layout" );
        FILEVARS = joinpath( PATHInfo, "variablesBRW.txt");
        writedlm( FILEVARS, TEXTVars );
        FILEPATHS = joinpath( PATHInfo, "Paths.jld" );
        PATHS = Dict(
            "PATHMain" => PATHMain,
            "PATHInfo" => PATHInfo,
            "PATHBRWs" => PATHBRWs
            )
        save( FILEPATHS, "PATHS", PATHS )
    else

        println( "You need another function for the: ", BRWversion, " version files, please consult the BRW manual" );
        Variables = Dict( );
        FILEPATHS = "";
    end

    return Variables, FILEPATHS

end
# -------------------------------------------------------------------------------------------------------- #
# using JLD
"""
SavePaths( FILEPATHS::String ) -> nothing
    Compares the name of the path variables in the workspace against the saved ones in the FILEPATHS file.
    If there are new paths, updates the file.
"""
function SavePaths( FILEPATHS::String )

    PATHS = load( FILEPATHS )[ "PATHS" ];
    vars = debug_list_vals( );
    pathsvars = String.( vars )[ findall( match.( r"PATH", String.( vars ) ) .!= nothing ) ];
    pathsfile = vcat( String.( keys( PATHS ) ), "FILEPATHS", "PATHS" );
    newpaths = pathsvars[ .!( pathsvars .∈ [ pathsfile ] ) ];

    if !isempty( newpaths )
        newdict = Dict( string( i ) => eval( Symbol( "$i" ) ) for i in newpaths );
        PATHS = merge( PATHS, newdict );
        save( FILEPATHS, "PATHS", PATHS );
        println( "Added ", keys( newdict ) );
    else
        println( "there is nothing new to add" )
    end

end
# -------------------------------------------------------------------------------------------------------- #
# using HDF5
"""
    Get_chunks( Variables::Dict{ String, Any }, Output_Chunks::String ) -> nothing
        Cuts the brw dataset into several more manageable segments
"""
function Get_chunks( Variables::Dict{ String, Any }, Output_Chunks::String )

    NRecFrames = Variables[ "NRecFrames" ];
    SamplingRate = Variables[ "SamplingRate" ];

    nChs = length( Variables[ "Layout" ] ); Σ = Variables[ "Raw" ];

    Σ = h5open( Variables[ "BRW" ], "r" )[ Σ ];

    ε, ω = size_segments( NRecFrames, SamplingRate );

    # number of spaces that occupies the number of seconds of duration of the experiment
    n_char = length( string( Int( floor( NRecFrames / SamplingRate ) ) ) );
    # number of bins
    n = size( ε, 1 );
    # number of spaces occupied by the number of chunks
    n_bins = length( string( n ) );
    #
    for i = 1:n # number of B to cut ( 1->4096, 1->ω )
        BIN = Array{ UInt16 }( undef, nChs, ω ); # preallocate
        #= values corresponding to the specific BIN.
         Channel 1,1 has the frames 1, 4097, 8193...etc =#
        β = collect( ( ε[ i, 1 ] - 1 ):ε[ i, 2 ] );
        for j = 1:ω
            # take those frames out of the vector Σ, and put them in array in BIN
            BIN[ :, j ] = Σ[ ( β[ j ] * nChs ) + 1:( nChs * β[ j + 1 ] ) ];
        end
        ini = lpad( string( Int(floor( ( β[ 1 ] )/SamplingRate ) ), "s" ), n_char + 1, "0" );
        eni = lpad( string( Int( floor( ( β[ end ] + 1 ) / SamplingRate ) ), "s" ), n_char + 1, "0" );
        bin_time = string( ini, "-", eni );
        BINname = string( "BIN", lpad( i, n_bins, "0" ), "_", bin_time, ".jld" );
        BINname = joinpath( Output_Chunks, BINname );
        save( BINname, "data", BIN );
    end
    close( Σ )
end
# -------------------------------------------------------------------------------------------------------- #
"""
    size_segments( NRecFrames::Int64, SamplingRate::Float64 ) -> ε::Vector, ω::Inf64
        Determines the ideal size of each chunck

"""
function size_segments( NRecFrames::Int64, SamplingRate::Float64 )
    if isinteger( SamplingRate ) # thus, there are chunks of 1 second
        if isinteger( NRecFrames/SamplingRate )
            n = Int.( NRecFrames/SamplingRate );
        end
    elseif isinteger( NRecFrames/floor( SamplingRate ) )
        n = Int.( NRecFrames/floor( SamplingRate ) );
    elseif isinteger( NRecFrames/ceil( SamplingRate ) )
        n = Int.( NRecFrames/ceil( SamplingRate ) );
    else # if not, a merequetengue is made
        div_T = div_ab( NRecFrames );
        div_sec = div_T/SamplingRate;
        # range of number of chunks to be cut, normally high to work at ease
        hi = 4; lo = 2; # segundos
        if !isempty( div_T )
            # finds one of the n-frames dividers within the range
            selected_divs = div_T[ findall( hi .>= div_sec .>= lo ) ];
        end
        if !isempty( selected_divs ) # if there was, grab the first one
            n = Int( NRecFrames/selected_divs[ 1 ] );
        else # if there wasn't, a default one.
            n = 60;
        end
    end
    #
    println( " The are ", ( n ), " segments, of ", ( ( NRecFrames / SamplingRate ) / n ) ," seconds each." );
    ω = ceil( Int, ( NRecFrames / n ) ); # number of final frames (chunk size in frames)
    ε = Array{ Int64 }( undef, n, 2 ); # preallocate
    ε[ :, 1 ] = collect( 1:ω:NRecFrames ); # start and
    ε[ :, 2 ] = ε[ :, 1 ] .+ ω .- 1; # end in frames of each chunk (to cut)
    if !isinteger( NRecFrames / n )
        ε[ end, 2 ] = NRecFrames; # end in frames of each chunk (to cut)
    end
    return ε, ω
end
# -------------------------------------------------------------------------------------------------------- #

"""
    div_ab( n::Int64, lo::Int = 1, hi::Int = n ) -> σ::Vector{Int64}
        Divisors of the number n between the values "lo" and "hi", if they are not defined then it takes from 1 to n
"""
function div_ab( n::Int, lo::Int = 1, hi::Int = n )
    ρ = collect( 1:floor( Int, sqrt( n ) ) ) ; # the numbers one by one, from the square root
    σ1 = findall( n.%ρ .== 0 ); # square root divisors ( remainder = 0 )
    σ2 = Int.( ( n ) ./ ( σ1 ) ); # Take out the pairs (of 100, 2-50, 10-10, etc.)
    σ = sort( unique( vcat( σ1, σ2 ) ) ); # remove duplicates, concatenate, sort
    aux1 = @isdefined lo;
    aux2 = @isdefined hi;
    if aux1 && aux2
        rn = σ[ findall( hi .>= σ .>= lo ) ];
        if isempty( rn )
            println(" There is no divisors of $n between $lo and $hi" )
        else
            return rn
        end
    else
        return σ
    end
end
# -------------------------------------------------------------------------------------------------------- #
"""
    Digital2Analogue( Variables::Dict{ String, Any }, Matrix::UInt16 ) -> Matrix{Float64}
        Conversion from raw data extracted from the brw file to voltage values (μV) acording to the equation
        Voltage = (RawData + ADCCountsToMV)*MVOffset
"""
function Digital2Analogue( Variables::Dict{ String, Any }, DigitalValue::Matrix{UInt16} )

    MVOffset = Variables[ "MVOffset" ];
    ADCCountsToMV = Variables[ "ADCCountsToMV" ];

    AnalogValue = @. MVOffset + ( DigitalValue * ADCCountsToMV )

    return AnalogValue

end
# -------------------------------------------------------------------------------------------------------- #
"""
    Digital2Analogue( Variables::Dict{ String, Any }, DigitalValue::UInt16 ) -> Float64
        Conversion from raw data extracted from the brw file to voltage values (μV) acording to the equation
        Voltage = (RawData + ADCCountsToMV)*MVOffset
"""
function Digital2Analogue( Variables::Dict{ String, Any }, DigitalValue::UInt16 )
    MVOffset = Variables[ "MVOffset" ]; ADCCountsToMV = Variables[ "ADCCountsToMV" ];
    AnalogValue = MVOffset + ( DigitalValue * ADCCountsToMV );
    return AnalogValue
end
# using Statistics
# using Distributions
# using BinningAnalysis
# -------------------------------------------------------------------------------------------------------- #
"""
    Noise-adaptive Optimal Thresholding
"""
@inline donoho( x ) =  ( median( abs.( x ) ) / 0.6745 );
# -------------------------------------------------------------------------------------------------------- #
@inline Sparsity( count_non_zeros::Int64, total_elements_of_A::Int64 ) = 1 - count_non_zeros/total_elements_of_A;
# -------------------------------------------------------------------------------------------------------- #
@inline Density( W::Vector{Any} ) = length( findall( ( length.( W ) ./ length( W ) ) .>= ( mean( length.( W )./length( W ) ) + 2*std( length.( W )./length( W ) ) ) ) );
# -------------------------------------------------------------------------------------------------------- #
"""
    Get_Groups( W::Vector{Int64} ) -> grupos::Vector{Any}, loose::( Vector{Any} )
       Groups of adjacent channels are formed from the initial indexes W. Those channels separated from the bulk of the others are considered loose ones.
"""
function Get_Groups( W::Vector{Int64} )
    grupos = [ ];
    for i in W
        _, vecinos = neighborgs( i, 1 );
        grupo = sort( intersect( vecinos, W ) );
        if isempty( grupo )
            push!( grupos, i )
        else
            push!( grupos, vcat( grupo, i ) )
        end
    end
    loose = grupos[ findall( length.( grupos ) .== 1 ) ];
    deleteat!( grupos, findall( length.( grupos ) .== 1 ) );
    a = length( grupos ); grupos = reverberation( grupos ); b = length( grupos );
    while a != b
        a = length( grupos ); grupos = reverberation( grupos ); b = length( grupos );
    end
    return grupos, loose
end
# -------------------------------------------------------------------------------------------------------- #
"""
    reverberation( grupos::Vector{Any} ) -> more_groups::Vector{Any}
        Intermediate step for Get_Groups function
"""
function reverberation( grupos::Vector{Any} )
    adjoining_channels = 0
    for i in grupos
        adjoining_channels = vcat( adjoining_channels, i )
    end
    adjoining_channels = adjoining_channels[ adjoining_channels .!= 0 ];
    adjoining_channels = unique( adjoining_channels );
    more_groups = [ ]
    for i in adjoining_channels
        temporal = 0
        for j = 1:length( grupos )
            if !isempty( findall( grupos[ j ] .== i ) )
                temporal = vcat( temporal, grupos[ j ] );
            end
        end
        temporal = temporal[ temporal .!= 0 ];
        new_group = sort( unique( temporal ) );

        if !isempty( new_group )
            if !isempty( more_groups )
                temp = last( more_groups );
                if !isequal( temp, new_group )
                    push!( more_groups, new_group );
                end
            else
                push!( more_groups, new_group );
            end
        end
    end
    return more_groups
end
# -------------------------------------------------------------------------------------------------------- #
"""
    neighborgs( C::Int64, d::Int64 ) ->
        -> A = Array( ( d*2 ) + 1, ( d * 2 ) + 1 ), v = vec( 2*( ( d * 2 ) + 1 ) - 1 );
        The d-neighborhood is calculated from the channel (C) as a center
        A = array where C is the center and is in chip order
        v = same neighboring channels as A but in vector form and without C ( 8 channels )
"""
function neighborgs( C::Int64, d::Int64 )
    Layout = reverse( reshape( collect( 1:4096 ), 64, 64 )', dims = 1 );
    x_c = findall( Layout .== C )[ ][ 2 ]; y_c = findall( Layout .== C )[ ][ 1 ];
    aux = [ ( x_c - d ),( x_c + d ), ( y_c - d ), ( y_c + d ) ]
    aux[ aux .< 1 ] .= 1; aux[ aux .> 64 ] .= 64;
    A = Layout[ aux[ 3 ]:aux[ 4 ], aux[ 1 ]:aux[ 2 ] ];
    v = vec( A )[ vec( A ) .!= C ];
    return A, v
end
# -------------------------------------------------------------------------------------------------------- #
"""
    MeanΔxCI( W::Vector, percent::Float64 ) -> xmean, Δx, C1, C2 (::Float64)
    Assumes a Normal distribution. Obtains the confidence interval with "percent" quantile,
    Returns mean, standard error, CI superior and CI inferior
"""
function MeanΔxCI( W::Vector, percent::Float64 )
    xmean, Δx = BinningAnalysis.jackknife( identity, W );
    d = Distributions.Normal( xmean, std( W ) );
    C1 = xmean + ( Δx * Distributions.quantile( d, percent ) );
    C2 = xmean - ( Δx * Distributions.quantile( d, percent ) );
    return xmean, Δx, C1, C2
end

# -------------------------------------------------------------------------------------------------------- #
# using Statistics
# using BinningAnalysis
"""
"""
function TimesPotentialDischarge( variability::Matrix{Float64} )
    R = copy( variability ); fill!( R, 0 );
    WT = zeros( size( variability, 2 ) );
    for j = 1:size( variability, 2 );
        W = variability[ :, j ];
        Todos, t = Thresholding( W );
        R[ :, j ] = sum( Todos, dims = 2 );
        xmean, xerror = jackknife( identity, R[ :, j ] );
        W = R[ :, j ];
        aux = union( findall( W .> xmean + std( W ) + xerror ), findall( W .< xmean - std( W ) - xerror ) );
        WT[ j ] = length( aux );
    end

    xmean, xerror = jackknife( identity, WT );
    Test04 = union( findall( WT .> xmean + 2*std( WT ) ), findall( WT .< xmean - 2*std( WT )  ) );
    return Test04
end

# -------------------------------------------------------------------------------------------------------- #
"""
"""
function SaturacionPositiva( Variables::Dict{String, Any}, data )

    nChs = length( Variables[ "Layout" ] );
    SamplingRate = Variables[ "SamplingRate" ];
    MinVolt = Float64( Variables[ "MinVolt" ] + 1 );
    channels = collect( 1:nChs );

    MinVoltSat = sum( Int.( data .== MinVolt ), dims = 2 );

    if length( findall( vec( MinVoltSat ) .== 0 ) ) == ( nChs - 1 )
        earth = setdiff( channels, findall( vec( MinVoltSat ) .== 0 ) )[ ];
        MinVoltSatChannels = [ ];
        data[ earth, : ] .= 0;
    else
        SatFramesLimit = ceil( Int, SamplingRate*0.005 );
        NineNinePercent = ceil( Int, binsize*0.99 );
        MinVoltSatChannels = findall( vec( NineNinePercent .>= MinVoltSat .>= SatFramesLimit ) );
        aux = findall( vec( MinVoltSat ) .>= NineNinePercent );
        if length( aux ) == 1
            earth = aux[ 1 ];
            data[ earth, : ] .= 0;
        else
            error( "there is something fishy here" );
        end
    end;
    return data, earth, MinVoltSatChannels
end

# -------------------------------------------------------------------------------------------------------- #
# using DSP
"""
"""
function FiltroMUAremez( Variables::Dict{String, Any}, canal::Vector{Float64} )

    SamplingRate = Variables[ "SamplingRate" ];
    lF = 300; fac = 10; HF = 3000;
    NYQ = floor( Int, SamplingRate / 2 );
    order = Int( floor( ( SamplingRate / lF ) / 5 ) );

    bpass = remez(
        ( order + 1 ), [ ( 0, lF - fac ) => 0, ( lF, HF ) => 1, ( HF + fac, NYQ ) => 0 ],
            Hz = SamplingRate );
    MUA = filtfilt( DSP.Filters.PolynomialRatio( bpass, [ 1.0 ] ), canal );
    return MUA
end

# -------------------------------------------------------------------------------------------------------- #
"""
"""
function ms2frames( time::Real, SamplingRate::Real )
    if time != 0
        x = ceil( Int, ( time * SamplingRate ) / 1000 );
    else
        x = 1
    end
    return x
end

# -------------------------------------------------------------------------------------------------------- #
# using SignalAnalysis
"""
"""
function SegmentsComplet( Temp02, antes, SegmentDuration, SamplingRate, limit, MUA )
    WindowIni = Temp02 .- ms2frames( antes, SamplingRate );
    total = antes*2 + SegmentDuration;
    WindowEnd = WindowIni .+ ms2frames( total, SamplingRate );
    SegmentDuration = ms2frames( SegmentDuration, SamplingRate );
    SegmentsCompletos = [ ];
    for aux = 1:length( Temp02 )
        wi = WindowIni[ aux ]; we = WindowEnd[ aux ];
        results = [ ]; wis = [ ];
        while ( we <= limit ) && wi >= 1 && ( wi + SegmentDuration < we)
            wbin = MUA[ wi:wi + SegmentDuration ];
            push!( results, energy( wbin ) );
            push!( wis, wi );
            wi = wi + 1;
        end
        if !isempty( results )
            wfinal = wis[ results .== maximum( results ) ][ 1 ];
            wfinal = [ wfinal, wfinal + SegmentDuration ];
            push!( SegmentsCompletos, wfinal )
        end
    end
    return Matrix( hcat( SegmentsCompletos... )' );
end

# -------------------------------------------------------------------------------------------------------- #
"""
"""
function SaturacionNegativaTemporal( Variables::Dict{String, Any}, data::VecOrMat{Float64} )

    SamplingRate = Variables[ "SamplingRate" ];
    MaxVolt = Float64( Variables[ "MaxVolt" ] - 1 );
    binsize = size( data, 2 );

    MaxVoltSatChannels = findall( vec( sum( Int.( data .== MaxVolt ), dims = 2 ) .!= 0 ) );
    MaxVoltFrames = [ ];
    for k = 1:length( MaxVoltSatChannels )
        push!( MaxVoltFrames, findall( data[ MaxVoltSatChannels[ k ], : ] .== MaxVolt ) );
    end

    maximo_permitido = ceil( Int, SamplingRate*0.3 );
    DiscartedChannels = MaxVoltSatChannels[ length.( MaxVoltFrames) .> maximo_permitido ];
    MaxSatReparables = setdiff( MaxVoltSatChannels, DiscartedChannels );
    dataTemp = copy( data );

    for k = 1:length( MaxSatReparables );

        Channel4Repair = MaxSatReparables[ k ];
        Frames4Repair = MaxVoltFrames[ MaxVoltSatChannels .==  Channel4Repair ][ ];
        FrameIndex = collect( 1:binsize );
        FramePool = setdiff( FrameIndex, Frames4Repair );
        FrameSample = sample( FramePool, length( Frames4Repair ) );
        VoltageSample = data[ Channel4Repair, FrameSample ];
        dataTemp[ Channel4Repair, Frames4Repair ] .= VoltageSample;

    end
    return dataTemp, DiscartedChannels

end

# -------------------------------------------------------------------------------------------------------- #
"""
"""
function ΔV( Variables::Dict{String, Any}, BIN::Matrix{Float64}, ΔT::Int64, descartados::Vector{Int64} )

    SamplingRate = Variables[ "SamplingRate" ][ 1 ];
    ΔT = ms2frames( ΔT, SamplingRate );

    STD = vec( std( ( BIN - circshift( BIN, ( 0, ΔT )  ) ), dims = 2 ) );
    STD2 = copy( STD ); STD2[ descartados ] .= 0;

    return STD, STD2

end

# -------------------------------------------------------------------------------------------------------- #
"""
"""
function Z0( X )

    Z = zeros( Int, 4096 );
    Z[ X ] = Z[ X ] .+ 1;
    Z = reverse( reshape( Z, 64, 64 )', dims = 1 );

    return Z
end

# -------------------------------------------------------------------------------------------------------- #
"""
"""
function ZW( X )
    Z = reverse( reshape( X, 64, 64 )', dims = 1 );
    return Z
end

# -------------------------------------------------------------------------------------------------------- #
"""
"""
function Zplot( Z, which, cm = :greys )
    if which == "0"
        Z = Z0( Z );
    elseif which == "W"
        Z = ZW( Z );
    end
    F = heatmap(
            Z,
            aspect_ratio = 1,
            c = cm,
            axis = ( [ ], false ),
            wsize = ( 400, 400 ),
            #cbar = :none
        );
    return F

end

# -------------------------------------------------------------------------------------------------------- #
# using HistogramThresholding
# using ImageContrastAdjustment
# using Suppressor
# using StatsBase
function Thresholding( W::VecOrMat{Float64} )
    W = copy( vec( W ) ); n = length( W );
    edges, conteo = HistogramThresholding.build_histogram( W, length( keys( countmap( W ) ) ) );

    t = zeros( 9 ); Todos = zeros( n, length( t ) );

    @suppress begin

        thr = find_threshold( UnimodalRosin( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 1 ] .= 1;
        t[ 1 ] = thr;
        thr = find_threshold( MinimumIntermodes( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 2 ] .= 1;
        t[ 2 ] = thr;
        thr = find_threshold( Intermodes( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 3 ] .= 1;
        t[ 3 ] = thr;
        thr = find_threshold( MinimumError( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 4 ] .= 1;
        t[ 4 ] = thr;
        thr = find_threshold( Moments( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 5 ] .= 1;
        t[ 5 ] = thr;
        thr = find_threshold( Otsu( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 6 ] .= 1;
        t[ 6 ] = thr;
        thr = find_threshold( Entropy( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 7 ] .= 1;
        t[ 7 ] = thr;
        thr = find_threshold( Balanced( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 8 ] .= 1;
        t[ 8 ] = thr;
        thr = find_threshold( Yen( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 9 ] .= 1;
        t[ 9 ] = thr;
    end

    return Int.( Todos ), t

end

# -------------------------------------------------------------------------------------------------------- #
# using Plots
function NuevePlots( Todos::Matrix{Int64} )

    AbstractImageBinarizationAlgorithm = [
        "UnimodalRosin",
        "MinimumIntermodes",
        "Intermodes",
        "MinimumError",
        "Moments",
        "Otsu",
        "Entropy",
        "Balanced",
        "Yen"
    ];

    m = 1
    xmessage = string( AbstractImageBinarizationAlgorithm[ m ], ": ", length( findall( Todos[ :, m ] .== 1 ) )," channels" );
    T1 = Zplot( Todos[ :, m ], "W" );
    T1 = plot!( title = xmessage, cbar = :none );

    m = 2
    xmessage = string( AbstractImageBinarizationAlgorithm[ m ], ": ", length( findall( Todos[ :, m ] .== 1 ) )," channels" );
    T2 = Zplot( Todos[ :, m ], "W" );
    T2 = plot!( title = xmessage, cbar = :none );

    m = 3
    xmessage = string( AbstractImageBinarizationAlgorithm[ m ], ": ", length( findall( Todos[ :, m ] .== 1 ) )," channels" );
    T3 = Zplot( Todos[ :, m ], "W" );
    T3 = plot!( title = xmessage, cbar = :none );

    m = 4
    xmessage = string( AbstractImageBinarizationAlgorithm[ m ], ": ", length( findall( Todos[ :, m ] .== 1 ) )," channels" );
    T4 = Zplot( Todos[ :, m ], "W" );
    T4 = plot!( title = xmessage, cbar = :none );

    m = 5
    xmessage = string( AbstractImageBinarizationAlgorithm[ m ], ": ", length( findall( Todos[ :, m ] .== 1 ) )," channels" );
    T5 = Zplot( Todos[ :, m ], "W" );
    T5 = plot!( title = xmessage, cbar = :none );

    m = 6
    xmessage = string( AbstractImageBinarizationAlgorithm[ m ], ": ", length( findall( Todos[ :, m ] .== 1 ) )," channels" );
    T6 = Zplot( Todos[ :, m ], "W" );
    T6 = plot!( title = xmessage, cbar = :none );

    m = 7
    xmessage = string( AbstractImageBinarizationAlgorithm[ m ], ": ", length( findall( Todos[ :, m ] .== 1 ) )," channels" );
    T7 = Zplot( Todos[ :, m ], "W" );
    T7 = plot!( title = xmessage, cbar = :none );

    m = 8
    xmessage = string( AbstractImageBinarizationAlgorithm[ m ], ": ", length( findall( Todos[ :, m ] .== 1 ) )," channels" );
    T8 = Zplot( Todos[ :, m ], "W" );
    T8 = plot!( title = xmessage, cbar = :none );

    m = 9
    xmessage = string( AbstractImageBinarizationAlgorithm[ m ], ": ", length( findall( Todos[ :, m ] .== 1 ) )," channels" );
    T9 = Zplot( Todos[ :, m ], "W" );
    T9 = plot!( title = xmessage, cbar = :none );

    TF = plot( T1, T2, T3, T4, T5, T6, T7, T8, T9, layout = ( 3, 3 ), wsize = ( 700, 700 ), titlefont = ( 8, "arial" ) );
    return TF
end

# -------------------------------------------------------------------------------------------------------- #
function FigureGroups( grupos::Vector, loose::Vector = [ ], cm = :twilight )
    Z = zeros( Int, 4096 );

    for i = 1:size( grupos, 1 )
        Z[ grupos[ i ] ] .= floor( Int, log( length( grupos[ i ] ) ) ) + 2 ;
    end
    Z[ Int.( loose ) ] .= 1
    Z = reverse( reshape( Z, 64, 64 )', dims = 1 )
    F = heatmap(
        Z,
        aspect_ratio = 1,
        c = cm,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        cbar = :none );
    return F
end

# -------------------------------------------------------------------------------------------------------- #
# using JuliaTools.neighborgs
function FillingHolesCrux( SatChannelsVar::Vector{Int64} )

    Z = reverse( reshape( 1:4096, 64, 64 )', dims = 1 );
    arriba = Z[ 64, 2:( end - 1 ) ]; abajo = Z[ 1, 2:( end - 1 ) ];
    izquierda = Z[ 2:( end - 1 ), 1 ]; derecha = Z[ 2:( end - 1 ), 64 ];
    esquinas = vcat( Z[ 1, 1 ], Z[ 64, 64 ], Z[ 64, 1 ], Z[ 1, 64 ] );
    bordes = vcat( arriba, abajo, izquierda, derecha, esquinas );
    bordes_noesquinas = vcat( arriba, abajo, izquierda, derecha );

    P = true
    while P
        X1 = setdiff( setdiff( setdiff( 1:4096, SatChannelsVar ) ), bordes );
        p1 = [ ];
        for i = 1:length( X1 )

            x = X1[ i ];
            A, _ = neighborgs( x, 1 );
            CruzH = [ A[ 2, 1 ], A[ 2, 3 ] ];
            CruzV = [ A[ 1, 2 ], A[ 3, 2 ] ];

            x1 = length( CruzH[ CruzH .∈ [ SatChannelsVar ] ] ) == 2;
            x2 = length( CruzV[ CruzV .∈ [ SatChannelsVar ] ] ) == 2;

            if ( x1 || x2 )
                push!( p1, x )
            end
        end

        X1 = intersect( arriba, setdiff( setdiff( 1:4096, SatChannelsVar ) ) );

        p2 = [ ];
        for i = 1:length( X1 )

            x = X1[ i ];
            A, _ = neighborgs( x, 1 );
                Cruz = [ A[ 2, 1 ], A[ 2, 3 ] ];
                x1 = length( Cruz[ Cruz .∈ [ SatChannelsVar ] ] ) == 2;
            if x1
                push!( p2, x )
            end
        end

        X1 = intersect( abajo, setdiff( setdiff( 1:4096, SatChannelsVar ) ) );

        p3 = [ ];
        for i = 1:length( X1 )

            x = X1[ i ];
            A, _ = neighborgs( x, 1 );
                Cruz = [ A[ 1, 1 ], A[ 1, 3 ] ];
                x1 = length( Cruz[ Cruz .∈ [ SatChannelsVar ] ] ) == 2;
            if x1
                push!( p3, x )
            end
        end

        X1 = intersect( vcat( izquierda, derecha ), setdiff( setdiff( 1:4096, SatChannelsVar ) ) );

        p4 = [ ];
        for i = 1:length( X1 )

            x = X1[ i ];
            A, _ = neighborgs( x, 1 );
                Cruz = [ A[ 1, 2 ], A[ 3, 2 ] ];
                x1 = length( Cruz[ Cruz .∈ [ SatChannelsVar ] ] ) == 2;
            if x1
                push!( p4, x )
            end
        end

         X1 = intersect( esquinas, setdiff( setdiff( 1:4096, SatChannelsVar ) ) );
         p5 = [ ];
         for i = 1:length( X1 )
             x = X1[ i ];
             _, v = neighborgs( x, 1 );

             if length( v[ v .∈ [ SatChannelsVar ] ] ) >= 2
                 push!( p5, x )
             end
         end
        posibles = [ ];
        posibles = vcat( Int.( p1 ), Int.( p2 ), Int.( p3 ), Int.( p4 ), Int.( p5 ) );

        if isempty( posibles )
            P = false
        else
            SatChannelsVar = vcat( posibles, SatChannelsVar );
        end

    end
    return Int.( SatChannelsVar )
end

end
