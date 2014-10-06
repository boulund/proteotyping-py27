#!/usr/bin/env julia
# Author: Fredrik Boulund
# Purpose: Identify bacterial species from short amino acid sequnces
# Date: 2014-07-24

using ArgParse
import PSMS

function parse_commandline()
    s = ArgParseSettings()
    s.description = "Bacterial proteotyping pipeline. Author: Fredrik Boulund"
    s.version = "0.0.1"
    s.add_version = true

    @add_arg_table s begin
        "FILE"
            help = "Input aa sequences."
            required = true
        "REF"
            help = "Reference database (FASTA) to compare against."
            required = false
        "--opt1"
            help = "An option with argument"
        "--debug"
            help = "A enable debug printouts"
            action = :store_true
    end

    add_arg_group(s, "output options")
    @add_arg_table s begin
        "--output", "-o"
            help = "Where to put the output"
            arg_type = String
            metavar = "OUTFILE"
            default = "output.txt"
    end

    return parse_args(s)
end


function runBLAT()
    database = "/collaborator/TTT/small_database_fixed_20140226.fasta"
    query = "test.fasta"
    output = "output.psl"
    run(`blat $database $query -prot $output`)
end




function main()
    parsed_args = parse_commandline()
    if parsed_args["debug"]
        println("Parsed args:")
        for (arg,val) in parsed_args
            println("  $arg:  $val")
        end
    end

    df = PSMS.parsePSMS(parsed_args["FILE"])
    PSMS.psmsToFasta(df, "test.fasta", remove_duplicates=true)

    runBLAT()

    
end


main()
