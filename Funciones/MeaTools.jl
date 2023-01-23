# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ #

module MeaTools

# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ #

using HDF5
using DelimitedFiles
using JLD
using Dates
using StatsBase
using DSP
using Statistics
using BinningAnalysis
using Plots
using MultivariateStats


# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ #

export GetVarsHDF5
export ChunkSizeSpace
export OneSegment
export Donoho
export EventsXBin
export EventsXChannel
export EventsXThrs
export STDΔV
export Ms2Frs
export MatrixFilter
export MUAremez
export Digital2Analogue
export DeSatMinMax
export AllMetrics
export WGeneral
export Zplot
export PosiblesDescargas

# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ #

# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    GetGroupsHDF5( BRW::HDF5.File, g::String ) -> GroupsN::Vector{String}
        Extract the groups form a BRW open file
"""
function GetGroupsHDF5( BRW::HDF5.File, g::String )
    GroupsN = [ ];
    try
        GroupsN = string.( g, "/", keys( BRW[ g ] ) );
    catch e
    end
    return GroupsN
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    GetAttrHDF5( BRW::HDF5.File, g::String ) -> AttrN::Vector{String}
        Extract the attributes form a BRW open file
        using HDF5
"""
function GetAttrHDF5( BRW::HDF5.File, g::String )
    AttrN = [ ];
    aux = attributes( BRW[ g ] );
    try
        AttrN = string.( g, "/", keys( aux ) );
    catch e
    end
    return AttrN
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    ExperimentSettings2Dict( Variables::Dict{Any, Any} ) -> Variables::Dict{Any, Any}
        Extracting the contents of the ExperimentSettings dictionary
"""
function ExperimentSettings2Dict( Variables::Dict{Any, Any} )
    ExperimentSettings = Variables[ "ExperimentSettings" ];
    t = split( ExperimentSettings, "\r\n" );
    t = replace.( t, "  " => "", "{" => "", "}" => "", '"' => "" );
    x = [ ];
    for i = 1:length( t )
        if !isempty( collect( eachmatch( r"[a-z]", t[ i ] ) ) )
            push!( x, true );
        else
            push!( x, false );
        end
    end
    t = t[ Bool.( x ) ]; t = split.( t, ": " );
    D = Dict( );
    for i in t
        if !( i[ 2 ] == "" )
            aux = i[ 2 ];
            try
                aux = replace( i[ 2 ], "," => " ", "[" => "", "[]" => "","]" => "", " " => "" )
            catch e
            end
            if ( aux != "" ) && ( aux != " " )
                aux = aux;
                try
                    aux = parse( Float64, aux );
                catch
                    aux = replace( aux, " " => "" );
                end
                D[ i[ 1 ] ] = aux;
            end
        end
    end
    delete!( Variables, "ExperimentSettings" );
    Variables = merge( Variables, D );
    return Variables
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    ExperimentDate2String( Variables::Dict{Any, Any} ) -> Variables::Dict{Any, Any}
        Extracting the date of the BRW creation
        using Dates
"""
function ExperimentDate2String( Variables::Dict{Any, Any} )
    X = Variables[ "ExperimentDateTime" ];
    Dt = split( X, ":" );
    Dt[ end ] = string( round( Int, parse( Float64, replace( Dt[ end ], r"[A-Z]" => "" ) ) ) );
    newDt = String( "" );
    for i in Dt
        newDt = string( newDt, ":", i );
    end
    newDt = newDt[ 2:end ];
    X = Dates.DateTime( newDt );
    Variables[ "ExperimentDateTime" ] = string( X );
    T = Dates.format( X, RFC1123Format );
    println( "Creation Date: ", T )
    return Variables
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    FindKeyRegex( key::String, D::Dict{Any, Any} ) -> Vector{String}
        Searchs the entries of the Dictionary D for the key word using Word Match.
"""
function FindKeyRegex( key::String, D::Dict{Any, Any} )
    key = Regex( ".+$key.+" );
    ok = keys( D );
    aux = match.( key, ok );
    aux01 = aux[ aux .!= nothing ];
    okok = [ ];
    if !isempty( aux01 )
        for i in aux01
            push!( okok, i.match );
        end
        return okok
    else
        println( "There is no match with that key word into the Dicctionary" );
    end
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
GetVarsHDF5( FILEBRW::String ) -> D::Dict{Any, Any}, FILEPATHS::String
        Extracts all contents from a HDF5 file of all versions, exept the Raw dataset, which only shows the string with the location
        using HDF5, DelimitedFiles, JLD, FindKeyRegex, ExperimentSettings2Dict, ExperimentDate2String, GetAttrHDF5, GetGroupsHDF5
"""
function GetVarsHDF5( FILEBRW::String )
    BRW = h5open( FILEBRW, "r" );
    Groups = keys( BRW );
    AllGroups = [ ];
    AllAtributes = keys( attributes( BRW ) );
    while !isempty( Groups )
        GROUPS = [ ];
        ATTR = [ ];
        for g in Groups
            if typeof( BRW[ g ] ) == HDF5.Dataset
                push!( AllGroups, g )
            else
                push!( GROUPS, GetGroupsHDF5( BRW, g ) );
                AllAtributes = vcat( AllAtributes, GetAttrHDF5( BRW, g ) );
            end
        end
        Groups = vcat( GROUPS... );
        push!( AllGroups, Groups );
    end
    AllGroups = vcat( AllGroups... );
    Types = Vector{ String }( undef, length( AllGroups ) );
    for g = 1:length( AllGroups )
        Types[ g ] = string( typeof( BRW[ AllGroups[ g ] ] ) );
    end
    AllDataSets = AllGroups[ Types .== "HDF5.Dataset" ];
    aux = zeros( Int, length( AllDataSets ) );
    for i = 1:length( AllDataSets )
        aux[ i ] = length( BRW[ AllDataSets[ i ] ] );
    end
    Raw = AllDataSets[ aux .== maximum( aux ) ][ 1 ];
    NoRaw = AllDataSets[ aux .!= maximum( aux ) ];
    aux = aux[ aux .!= maximum( aux ) ];
    NoRaw = NoRaw[ aux .!= 0 ];
    D = Dict( );
    e = [ ];
    for g in NoRaw
        try
            D[ g ] = Float64.( read( BRW[ g ] ) );
        catch e
            D[ g ] = read( BRW[ g ] );
        end
    end
    for g in keys( D )
        if length( D[ g ] ) .== 1
            D[ g ] = D[ g ][ 1 ];
        end
        try
            D[ g ] = Int64( D[ g ] );
        catch e
        end
    end
    for a in AllAtributes
        aux00 = basename( a );
        aux01 = dirname( a );
        if isempty( aux01 )
            D[ a ] = read_attribute( BRW, aux00 );
        else
            D[ a ] = read_attribute( BRW[ aux01 ], aux00 );
        end
    end
    D[ "Raw" ] = Raw;
    BRWsize = ( ( stat( FILEBRW ).size ) / 1000000 ) / 1024;
    D[ "BRWname" ] = BRW.filename;
    D[ "BRWsizeGB" ] = BRWsize;
    aux0 = FindKeyRegex( "Ch", D );
    aux1 = [ ];
    for i in aux0
        try
            push!( aux1, size( D[ i ], 1 ) );
        catch e
            push!( aux1, 0 );
        end
    end
    nChs = size( D[ aux0[ aux1 .== maximum( aux1 ) ][ ] ], 1 );
    D[ "nChs" ] = nChs;
    aux0 = FindKeyRegex( "Std", D );
    aux1 = [ ];
    for i in aux0
        try
            push!( aux1, size( D[ i ], 1 ) );
        catch e
            push!( aux1, 0 );
        end
    end
    STD = D[ aux0[ aux1 .== maximum( aux1 ) ][ ] ];
    D[ "STD" ] = STD;
    NewVars = Dict( );
    K = keys( D );
    for k in K
        NewVars[ basename( k ) ] = D[ k ];
    end
    D = NewVars;
    if ( "ExperimentSettings" in keys( D ) )
        D = ExperimentSettings2Dict( D );
        D = ExperimentDate2String( D );
    end
    X = String.( keys( D ) )[ values( D ) .== "null" ];
    for i in X
        delete!( D, i );
    end
    X = [ ];
    for i in keys( D )
        try
            if isempty( D[ i ] )
                push!( X, i );
            end
        catch e
        end
    end
    for i in X
        delete!( D, i );
    end
    x = [ ];
    for i in values( D )
        if length( i ) <= 100
            push!( x, 1 )
        else
            push!( x, 0 )
        end
    end
    NK = string.( keys( D ) )[ Bool.( x ) ];
    TXT = Dict( ); for i in NK; TXT[ i ] = D[ i ]; end
    K = keys( D );
    NewVars = Dict( );
    for k in K
        NewVars[ basename( k ) ] = D[ k ];
    end
    D = NewVars;
    PATHMAIN = split( FILEBRW, "." )[ 1 ]; mkpath( PATHMAIN );
    PATHINFO = joinpath( PATHMAIN, "Info" ); mkpath( PATHINFO );
    FILEVARS = joinpath( PATHINFO, "Variables.jld" );
    FILEVARSTXT = joinpath( PATHINFO, "Variables.txt" );
    FILEPATHS = joinpath( PATHINFO, "Paths.jld" );
    PATHBRWs = dirname( FILEBRW );
    PATHS = Dict(
        "PATHMAIN" => PATHMAIN,
        "PATHINFO" => PATHINFO,
        "PATHBRWs" => PATHBRWs
        );
    writedlm( FILEVARSTXT, TXT );
    save( FILEVARS, "Variables", D );
    save( FILEPATHS, "PATHS", PATHS );
    println( "You are now working on the new main path: ", PATHMAIN );
    println( "With the file: ")
    print( basename( BRW.filename ), " : ", D[ "Description" ] );
    println( " HDF5 file size: $BRWsize GB" );
    cd( PATHMAIN )
    close( BRW )
    return D, FILEPATHS
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    ChunkSizeSpace( Variables::Dict{Any, Any}, limupper::Real ) -> σ::Real
    Using the maximum size ( limupper in GB ) given, determines the number of segmentes σ to extract from the entire dataset
"""
function ChunkSizeSpace( Variables::Dict{Any, Any}, limupper::Real )
    NRecFrames = Variables[ "NRecFrames" ];
    σ = 5
    finalsize = Variables[ "BRWsizeGB" ] / σ;
    while ( finalsize > limupper )
        σ = σ + 1
        finalsize = Variables[ "BRWsizeGB" ] / σ;
    end
    while NRecFrames % σ != 0
        σ = σ + 1
    end
    println( "$σ segments of $finalsize GB each one" )
    return σ
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    OneSegment( Variables::Dict{Any, Any}, n::Int64, nSegments::Int ) -> BIN::::Matrix{Float64}
        Takes just one segment of the BRW file
        using HDF5
"""
function OneSegment( Variables::Dict{Any, Any}, n::Int64, nSegments::Int )
    nChs = Variables[ "nChs" ];
    NRecFrames = Variables[ "NRecFrames" ];
    binlenght = Int( NRecFrames / nSegments );
    Raw = h5open( Variables[ "BRWname" ], "r" )[ Variables[ "Raw" ] ];
    BIN = Array{ UInt16 }( undef, nChs, binlenght );
    for frame = ( ( ( n - 1 ) * binlenght ) + 1 ): binlenght * n
        BIN[ :,
            ( frame - ( binlenght*( n - 1 ) ) ) ] .=
                Raw[ ( ( ( frame - 1 ) * nChs ) + 1 ): nChs * frame ];
    end
    return BIN
    close( Raw )
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    Digital2Analogue( Variables::Dict{ Any, Any }, DigitalValue::Matrix{UInt16} )
        -> BIN::Matrix{Float64}
        Voltage Convertion
"""
function Digital2Analogue( Variables::Dict{ Any, Any }, DigitalValue::Matrix{UInt16} )
    SignalInversion = Variables[ "SignalInversion" ];
    MinVolt = Variables[ "MinVolt" ];
    MaxVolt = Variables[ "MaxVolt" ];
    BitDepth = Variables[ "BitDepth" ];
    MVOffset = SignalInversion*MinVolt;
    ADCCountsToMV = ( SignalInversion * ( MaxVolt - MinVolt ) ) / ( 2^BitDepth );
    BIN = @. MVOffset + ( DigitalValue * ADCCountsToMV )
    return BIN
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
 DeSatMinMax( Variables::Dict{Any, Any}, BIN::Matrix{Float64}, ChPercent = 15, FrPercent = 20 )
    -> BIN::Matrix{Float64}, DiscardedFrames::Vector{Int64}, DiscardedChannels::Vector{Int64},
    FrsForRepair::Vector{Int64}, ChsForRepair::Vector{Int64}
    using StatsBase
"""
function DeSatMinMax( Variables::Dict{Any, Any}, BIN::Matrix{Float64}, ChPercent = 15, FrPercent = 20 )
    MaxVolt = Variables[ "MaxVolt" ] - 100;
    MinVolt = Variables[ "MinVolt" ] + 100;
    nChs = size( BIN, 1 );
    nFrs = size( BIN, 2 );
    LCh = ceil( Int, nFrs * ChPercent / 100 );
    LFr = ceil( Int, nChs * FrPercent / 100 );
    #
    ChFrMax = getindex.( findall( BIN .>= MaxVolt ), [ 1 2 ] );
    ChFrMin = getindex.( findall( BIN .<= MinVolt ), [ 1 2 ] );
    SatChs = zeros( Int, nChs, 3 );
    for i = 1:size( ChFrMax, 1 )
        SatChs[ ChFrMax[ i, 1 ], 1 ] = SatChs[ ChFrMax[ i, 1 ], 1 ] + 1;
    end
    for i = 1:size( ChFrMin, 1 )
        SatChs[ ChFrMin[ i, 1 ], 2 ] = SatChs[ ChFrMin[ i, 1 ], 2 ] + 1;
    end
    SatChs[ :, 3 ] .= SatChs[ :, 1 ] .+ SatChs[ :, 2 ];
    ChsForRepair = findall( SatChs[ :, 3 ] .!= 0 );
    DiscardedChannels = ChsForRepair[ SatChs[ ChsForRepair, 3 ] .>= LCh ];
    PossibleEarth = ChsForRepair[ SatChs[ ChsForRepair, 3 ] .>= nFrs*0.8 ];
    ChsForRepair = setdiff( ChsForRepair, DiscardedChannels );
    #
    SatFrs = zeros( Int, nFrs, 3 );
    for i = 1:size( ChFrMax, 1 )
        SatFrs[ ChFrMax[ i, 2 ], 1 ] = SatFrs[ ChFrMax[ i, 2 ], 1 ] + 1;
    end
    for i = 1:size( ChFrMin, 1 )
        SatFrs[ ChFrMin[ i, 2 ], 2 ] = SatFrs[ ChFrMin[ i, 2 ], 2 ] + 1;
    end
    SatFrs[ :, 3 ] .= SatFrs[ :, 1 ] .+ SatFrs[ :, 2 ];
    FrsForRepair = findall( SatFrs[ :, 3 ] .!= 0 );
    DiscardedFrames = FrsForRepair[ SatFrs[ FrsForRepair, 3 ] .>= LFr ];
    if isempty( DiscardedFrames )
            FrsForRepair = [ ];
    else
        FrsForRepair = setdiff( FrsForRepair, DiscardedFrames );
    end
    for ch in ChsForRepair
        Frames = vcat( findall( BIN[ ch, : ] .>= MaxVolt ), findall( BIN[ ch, : ] .<= MinVolt ) );
        Frames = setdiff( Frames, DiscardedFrames );
        BuenosFrames = setdiff( 1:nFrs, Frames );
        BIN[ ch, Frames ] = BIN[ ch, sample( BuenosFrames, length( Frames ) ) ];
    end
    return BIN, DiscardedFrames, DiscardedChannels, FrsForRepair, ChsForRepair, PossibleEarth
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    MUAremez( Variables::Dict{Any, Any}, channel::Vector{Float64} ) -> MUA::Vector{Float64}
    using DSP
"""
function MUAremez( Variables::Dict{Any, Any}, channel::Vector{Float64} )
    lF = 300;
    fac = 10;
    HF = 3000;
    SamplingRate = Variables[ "SamplingRate" ];
    NYQ = floor( Int, SamplingRate / 2 );
    order = Int( floor( ( SamplingRate / lF ) / 5 ) );
    bpass = remez(
        ( order + 1 ), [ ( 0, lF - fac ) => 0, ( lF, HF ) => 1, ( HF + fac, NYQ ) => 0 ],
            Hz = SamplingRate );
    MUA = filtfilt( DSP.Filters.PolynomialRatio( bpass, [ 1.0 ] ), channel );
    return MUA
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    MatrixFilter( Variables::Dict{Any, Any}, data::Matrix{Float64} ) -> datafilt::Matrix{Float64}
    using FiltroMUAremez
"""
function MatrixFilter( Variables::Dict{Any, Any}, data::Matrix{Float64} )
    datafilt = copy( data );
    for k = 1:size( data, 1 )
        channel = data[ k, : ];
        MUA = MUAremez( Variables, channel );
        datafilt[ k, : ] = MUA;
    end
    return datafilt
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    Ms2Frs( time::Real, Variables::Dict{Any, Any} ) -> x::Float64
"""
function Ms2Frs( time::Real, Variables::Dict{Any, Any} )
    SamplingRate = Variables[ "SamplingRate" ];
    if time != 0; x = ceil( Int, ( time * SamplingRate ) / 1000 ); else; x = 1; end
    return x
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
STDΔV( Variables::Dict{Any, Any}, BIN::Matrix{Float64}, ΔT::Real ) -> STD::Vector{Float64}
        using Ms2Frs
"""
function STDΔV( Variables::Dict{Any, Any}, BIN::Matrix{Float64}, ΔT::Real )
    ΔT = Ms2Frs( ΔT, Variables );
    STD = vec( std( ( BIN - circshift( BIN, ( 0, ΔT )  ) ), dims = 2 ) );
    return STD
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    EventsXThrs( channel::Vector{Float64}, thr::Real, parameters::Dict{String, Int64} )
    -> index::Vector{Int64}
    Index of detected suprathreshold events with a pre-established threshold
    using EventsXBin
"""
function EventsXThrs( channel::Vector{Float64}, thr::Real, parameters::Dict{String, Int64} )
    index = EventsXBin( channel, thr, parameters );
    return index
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    EventsXChannel( channel::Vector{Float64}, parameters::Dict{String, Int64} )
    -> thrs::Vector{Float64}, index_real::Vector{Inf64}
    using Donoho, EventsXBin
"""
function EventsXChannel( channel::Vector{Float64}, parameters::Dict{String, Int64} )
    window = parameters[ "window" ];
    bit = parameters[ "bit" ];
    σ = parameters[ "cte" ];
    index_real = [ ];
    i = 1;
    I = Int( ( ( ( i - 1 ) * bit ) + 1 ) ); J = Int( I + window - 1 );
    thrs = [ ];
    thr = 0;
    while J <= ( length( channel ) - Int( window - 1 ) )
        bin = channel[ I:J ];
        if !iszero( bin )
            thr = -1*σ*abs( Donoho( bin ) );
            index_parcial = EventsXBin( bin, thr, parameters );
            if !isempty( index_parcial )
                index_real = vcat( index_real, ( index_parcial .+ I .- 1 ) );
            end
        end
        i = i + 1;
        I = Int( ( ( ( i - 1 ) * bit ) + 1 ) ); J = Int( I + window - 1 );
    end
    index_real = unique( index_real );
    if !isempty( index_real )
        push!( thrs, thr )
    else
        push!( thrs, [0] )
    end
    return vcat( thrs... ), index_real
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    EventsXBin( bin::Vector{Float64}, thr::Real, parameters::Dict{String, Int64} ) -> Index::Vector{Int64}
"""
function EventsXBin( bin::Vector{Float64}, thr::Real, parameters::Dict{String, Int64} )
    distance = parameters[ "distance" ];
    ST = findall( bin .<= thr ); # eventos que pasan el umbral
    index_parcial = [ ];
    if !isempty( ST )
        a = 1;
        while a == 1
            distances = diff( ST ); # distancia entre ellos
            nears = findall( distances .<= distance ) .+ 1; # cuales estan cerca
            if isempty( nears )
                a = 0;
            else
                remove = zeros( Int, size( nears, 1 ) );
                for i = 1:size( nears, 1 )
                    # si el primero es menor que el segundo, quita el segundo
                    if isless( ( bin[ ST[ nears[ i ] ] ] ),
                            ( bin[ ST[ nears[ i ] - 1 ] ] ) )
                        remove[ i ] = nears[ i ] - 1;
                    else
                        remove[ i ] = nears[ i ];
                    end
                end
                if size( remove, 1 ) > 1
                    ST[ unique( remove ) ] .= 0;
                    filter!( x -> x != 0, ST );
                else
                    ST = ST[ Bool.( ST .!= ST[ remove[ 1 ] ] ) ];
                end
            end
        end
        push!( index_parcial, ST )
    else
        push!( index_parcial, [ ] )
    end
    return sort( unique( vcat( index_parcial... ) ) )
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    Donoho( x::Vector ) -> thr::Float64
        Noise-adaptive Optimal Thresholding
"""
Donoho( x ) = ( median( abs.( x ) ) / 0.6745 );
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
AllMetrics( Variables::Dict{Any, Any}, BIN::Matrix{Float64}, ΔT::Real, window::Real, bit::Real, distance::Real, cte::Real, thr::Real ) -> Ws::Matrix{Float64}, AllEvents_RS_RD_FS_FD::Vector{Vector{Any}}

"""
function AllMetrics( Variables::Dict{Any, Any}, BIN::Matrix{Float64}, ΔT::Real, window::Real, bit::Real, distance::Real, cte::Real, thr::Real )
    nChs = Variables[ "nChs" ];
    channels = 1:nChs;
    Ws = zeros( Float64, nChs, 9*2 );
    parameters = Dict(
        "window"    => Ms2Frs( window, Variables ),
        "bit"      => Ms2Frs( bit, Variables ),
        "distance"  => Ms2Frs( distance, Variables ),
        "cte"       => cte
    );
    FiltBIN = MatrixFilter( Variables, BIN );
    # 1 Desviacion estandar del cambio de voltaje con Δt variable
    Ws[ :, 1 ] = STDΔV( Variables, BIN, ΔT );
    Ws[ :, 2 ] = STDΔV( Variables, FiltBIN, ΔT );
    # 2 Numero de valores de voltaje presentes en cada canal
    [ Ws[ k , 3 ] = length( unique( round.( BIN[ k, : ], digits = 2 ) ) ) for k in channels ];
    [ Ws[ k , 4 ] = length(
        unique( round.( FiltBIN[ k, : ], digits = 2 ) ) ) for k in channels ];
    # 3 y 4
    EventsFiltDin = [ ]; EventsRawDin = [ ];
    for i in channels
        _, IndexFilt = EventsXChannel( FiltBIN[ i, : ], parameters );
        _, IndexRaw = EventsXChannel( BIN[ i, : ], parameters );
        push!( EventsFiltDin, IndexFilt );
        push!( EventsRawDin, IndexRaw );
    end
    EventsRawStat = [ ]; EventsFiltStat = [ ];
    for i in channels
        IndexRaw = EventsXThrs( BIN[ i, : ], thr, parameters );
        IndexFilt = EventsXThrs(  FiltBIN[ i, : ], thr, parameters );
        push!( EventsFiltStat, IndexFilt );
        push!( EventsRawStat, IndexRaw );
    end
    Ws[ :, 5 ] = length.( EventsRawStat );
    Ws[ :, 6 ] = length.( EventsRawDin );
    Ws[ :, 7 ] = length.( EventsFiltStat );
    Ws[ :, 8 ] = length.( EventsFiltDin );
    AllEvents_RS_RD_FS_FD = [ EventsRawStat, EventsRawDin, EventsFiltStat, EventsFiltDin ];
    # 5 Desviacion estandar directa
    Ws[ :, 9 ] = std( BIN, dims = 2 );
    Ws[ :, 10 ] = std( FiltBIN, dims = 2 );
    # 7. medians
    Ws[ :, 13 ] = median( BIN, dims = 2 );
    Ws[ :, 14 ] = median( FiltBIN, dims = 2 );
    # 8. Global threshold
    [ Ws[ k , 15 ] = Donoho( BIN[ k, : ] ) for k in channels ];
    [ Ws[ k , 16 ] = Donoho( FiltBIN[ k, : ] ) for k in channels ];
    # 6, means and 9, error
    for k in channels
        Ws[ k, 11 ], Ws[ k, 17 ] = jackknife( identity, BIN[ k, : ] );
        Ws[ k, 12 ], Ws[ k, 18 ] = jackknife( identity, FiltBIN[ k, : ] );
    end
    return Ws, AllEvents_RS_RD_FS_FD, Float16.( FiltBIN )
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
"""
    Z0( X::VecOrMat, nChs::Int64 ) -> Z::Matrix{Int64}
        using Plots
"""
function Z0( X::VecOrMat, nChs::Int64 )
    X = Int.( vec( X ) );
    Z = zeros( Int, nChs );
    n = Int( sqrt( nChs ) );
    Z[ X ] .= Z[ X ] .+ 1;
    Z = reverse( reshape( Z, n, n )', dims = 1 );
    return Int.( Z )
end
"""
    ZW( X::VecOrMat ) -> Z::Matrix{typeof(X)}
        using Plots
"""
function ZW( X::VecOrMat )
    X = vec( X );
    n = Int( sqrt( length( X ) ) );
    Z = reverse( reshape( X, n, n )', dims = 1 );
    return Z
end
"""
    Zplot( Z::Matrix, which::String, cm = :greys, nChs = 4096 ) -> F::Plot
        using Plots, Z0, ZW
"""
function Zplot( Z::VecOrMat, which::String, cm = :greys, nChs = 4096 )
    if which == "0"
        Z = Z0( Z, nChs );
    elseif which == "W"
        Z = ZW( Z );
    end
    F = heatmap( Z, aspect_ratio = 1, c = cm, axis = ( [ ], false ), wsize = ( 400, 400 ) );
    return F
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
function WGeneral( Wn::Matrix{Float64}, n::Int64, σ::Int64, mid::Float64 )
    titles = [
        "std ΔV, log, N, Raw", "std ΔV, log, N, Filt", "# Voltage, log, N, Raw",
        "# Voltage, log, N, Filt", "EventsStat, log, N, Raw", "EventsDin, log, N, Raw",
        "EventsStat, log, N, Filt", "EventsDin, log, N, Filt", "Std, log, N, Raw",
        "Std, log, N, Filt", "mean, log, N, Raw", "mean, log, N, Filt",
        "median, log, N, Raw", "median, log, N, Filt", "threshold, log, N, Raw",
        "threshold, log, N, Filt", "error, log, N, Raw", "error, log, N, Filt",
    ];
    indexes = collect( ( ( n - 1 )*σ + 1 ) : σ*n );
    W = vec( sum( Wn[ :, indexes ], dims = 2 ) );
    W = abs.( round.( W, digits = 3 ) ); W0 = round.( log.( W ), digits = 3 );
    W0N = copy( W0 );
    p = round( minimum( W0[ W0 .!= -Inf ] ), digits = 3 ); W0N[ W0N .== -Inf ] .= p;
    cm = reverse( colormap( "RdBu", length( countmap( W0N ) ); mid = mid, b = 0.01 ) );
    P = Zplot( W0N, "W", cm ); P = plot!( title = titles[ n ], titlefontsize = 8 );
    return P, W0N
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
function GetM( W, v2 )
    M = Array{Float64}( undef, v2, 8 );
    M[ :, 1 ] = median( W, dims = 1 );
    M[ :, 2 ] = mean( W, dims = 1 );
    M[ :, 3 ] = std( W, dims = 1 );
    M[ :, 4 ] = mean(W, dims = 1) .+ 2*std( W, dims = 1 );
    N = [ ]
    for i in 1:size( W, 2 )
        push!( N, Donoho( W[ :, i ] ) );
    end
    M[ :, 5 ] = N;
    N = [ ]
    for i in 1:size( W, 2 )
        push!( N, ( Donoho( W[ :, i ] ) + sqrt( 2*log( length( W[ :, i ] ) ) ) ) );
    end
    M[ :, 6 ] = N;
    M[ :, 7 ] = var( W, dims = 1 );
    N = [ ]
    for i in 1:size( W, 2 )
        push!( N, mode( W[ :, i ] ) );
    end
    M[ :, 8 ] = N;
    return M
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
function ForSomeReason( AllWs::Dict{String, Any}, nM )
    Wn = AllWs[ "Wn" ];
    W0NTotal = AllWs[ "W0NTotal" ];
    v2 = Int( size( Wn, 2 ) / size( W0NTotal, 2 ) );
    indexes = ( ( nM - 1 )*v2 + 1 ) : ( nM*v2 );
    W = Wn[ :, indexes ]; M = GetM( W, v2 );
    R = fit( ICA, M, 2 ); NR = R.W[ :, 1 ];
    _, xerror = jackknife( identity, NR );
    thr1 = median( NR ) - 2 * xerror;
    thr2 = median( NR ) + 2 * xerror;
    bins = sort( union( findall( NR .>= thr2 ), findall( NR .<= thr1 ) ) );
    return bins
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·• #
function PosiblesDescargas( AllWs, nM = 3 )
    test = [ ]
    for i in 1:10
        bins = ForSomeReason( AllWs, nM );
        push!( test, bins )
    end
    bines = unique( vcat( test... ) );
    return bines
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ #

end # ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
