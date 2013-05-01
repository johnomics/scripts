using GZip

import Base.length

type VCFReader
    v::IO
    VCFReader(filename::String)=new(gzopen(filename))
end

function start(vcf::VCFReader)
    readline(vcf.v)
    return 0
end

next(vcf::VCFReader,snps)= readline(vcf.v), snps+1

done(vcf::VCFReader,snps) = eof(vcf.v)

type Scaffold
    name::String
    length::Int
    snps::Vector{String}
    Scaffold(scf::String,len::Int)=new(scf,len,String[])
    Scaffold(scf::String)=Scaffold(scf,0)
    Scaffold()=Scaffold("")
end

type Genome
    vcf::VCFReader
    length::Int
    scaffolds::Array{Scaffold}
    partsz::Int

    function Genome(filename,np::Int)
        l=0
        s::Array{Scaffold}={}
        v=VCFReader(filename)
        for header in v
            if ismatch(r"^#CHROM",header)
                break
            end
            
            m=match(r"^##contig=<ID=(.+),length=(\d+)>$",header)
            if m != nothing
                scfname = m.captures[1]
                scflen = int(m.captures[2])
                push!(s,Scaffold(scfname,scflen))
                l += scflen
            end
        end
        new(v,l,s,int(l/np))
    end
end

function outputScaffold(s::Scaffold)
    println("$(s.name)\t$(s.length)")
end

function outputScaffolds(g::Genome)
    for s in g.scaffolds
        outputScaffold(s)
    end
end

function length(g::Genome)
    g.length
end