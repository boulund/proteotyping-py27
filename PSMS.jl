#!/usr/bin/env julia
# Author: Fredrik Boulund
# Purpose: Read PSMS output txt files
# Date: 2014-07-24

module PSMS

using DataFrames


function parsePSMS(filename::String)
    return readtable(filename, separator='\t')
end


function psmsToFasta(psms_data, filename; remove_duplicates=false)
    out = open(filename, "w")

    if remove_duplicates
        sequences = Set(map(uppercase, psms_data[4]))
    else
        sequences = psms_data[4]
    end

    counter = 1
    for seq in sequences
        println(out, ">", counter)
        println(out, seq)
        counter += 1
    end

end


end
