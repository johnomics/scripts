#!/usr/bin/env julia

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

function next(vcf::VCFReader, snps)
    return readline(vcf.v),snps
end

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
    filename = args[1]
    nprocs = int(args[2])
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