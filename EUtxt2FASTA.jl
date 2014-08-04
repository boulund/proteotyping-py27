#!/usr/bin/env julia
# Fredrik Boulund
# 2014-08-04
# Convert a psms output file to FASTA
import PSMS
function main()

    if length(ARGS)<2
        println("usage: PSMS.jl psms.txt outfile.fasta")
        return
    end

    df = PSMS.parsePSMS(ARGS[1])
    PSMS.psmsToFasta(df, ARGS[2], remove_duplicates=false)
end

main()
