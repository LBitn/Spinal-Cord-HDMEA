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

export div_ab
export size_segments
export Get_chunks
export SavePaths
export debug_list_vals
export VariablesBRW
export searchdir
export Digital2Analogue
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
# using DataFrames
"""
debug_list_vals( m::Module = Main )-> vs::Vector{Symbol}
    Creates a list of the names of the variables presents in the current workspace
"""
function debug_list_vals( m::Module = Main )
    vs = [ name for name in sort!( names( m ) ) if isdefined( m, name ) ];
    return vs
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

end
