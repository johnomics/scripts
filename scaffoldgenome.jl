#!/usr/bin/env julia

using ArgParse
using GZip
import Base.start, Base.done, Base.next



function main(args)
    s=ArgParseSettings("Scaffold genome sequence from genetic marker information\n")
    
    @add_arg_table s begin
        "--procs", "-p"
            arg_type = Int
            default = 1
            help = "Number of processors"
        "--maxsnps", "-m"
            arg_type = Int
            default = 0
            help = "Maximum SNPs to process"
        "filename"
            arg_type = String
            required = true
            help = "VCF file containing marker SNPs"
    end
    
    parsed_args = parse_args(args, s)
    filename = parsed_args["filename"]
    np = parsed_args["procs"]
    maxsnps = parsed_args["maxsnps"]
    
    addprocs(np-1)
    
    # Require has to be run after all processors are added
    require("/whale-data/jd626/scripts/vcf.jl")
    
    snps = 0
    
    # Process header
    g = Genome(filename,np)
    outputScaffolds(g)
    println(length(g))
end

#    curproc = 1
#    for snp in eachline(vcf)
#        if ismatch(r"^#", snp)
#            continue
#        end
#        r[curproc] = @spawnat curproc push!(lines,snp)
#        curproc+=1
#        if curproc > np curproc = 1 end
#        snps +=1
#        if maxsnps > 0 && maxsnps <= snps
#            break
#        end
#        if snps % 1000 == 0 print('.') end
#        if snps % 10000 == 0 print("$snps\n") end
#    end
#end

main(ARGS)