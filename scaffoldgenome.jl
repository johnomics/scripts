#!/usr/bin/env julia

using ArgParse
using GZip
import Base.start, Base.done, Base.next


type VCFReader
    v::IO
    VCFReader(filename::String)=new(gzopen(filename))
end

function start(vcf::VCFReader)
    readline(vcf.v)
    return 0
end

next(vcf::VCFReader, snps)=readline(vcf.v),snps

done(vcf::VCFReader,snps) = eof(vcf.v)

type Scaffold
    name::String
    snps::Vector{String}
    Scaffold(scf::String)=new(scf,String[])
    Scaffold()=Scaffold("")
end

function processScaffold(scf::Scaffold)
    return length(scf.snps)
end

function main(args)
    s=ArgParseSettings("Scaffold genome sequence from genetic marker information\n")
    
    @add_arg_table s begin
        "--procs", "-p"
            arg_type = Int
            default = 1
            help = "Number of processors"
        "filename"
            arg_type = String
            required = true
            help = "VCF file containing marker SNPs"

    end
    
    parsed_args = parse_args(args, s)
    filename = parsed_args["filename"]
    nprocs = parsed_args["procs"]

    addprocs(nprocs-1)
    vcf = VCFReader(filename)
    snps = 0
    scf = Scaffold()
    scfrefs = Dict{String,Int}()

    for (snp) in vcf
        if (ismatch(r"^#", snp)) continue end
        snps += 1
        snpscf = split(snp,"\t")[1]
        if (snpscf != scf.name )
            if (scf.name != "")
                scfrefs[scf.name] = processScaffold(scf)
            end
            scf=Scaffold(snpscf)
        end
        push!(scf.snps,snp)
        if (snps %  10000 == 0)
            print(".")
        end
        if (snps % 100000 == 0)
            println(snps)
        end
    end
    
    for (scf) in sort(keys(scfrefs))
        println("$scf\t$(scfrefs[scf])")
    end
end

main(ARGS)