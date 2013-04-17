#!/usr/bin/env julia

using GZip
import Base.start, Base.done, Base.next


type VCFReader
    v::IO
    function VCFReader(filename::String)
        new(gzopen(filename))
    end
end

function start(vcf::VCFReader)
    readline(vcf.v)
    return 0
end

function next(vcf::VCFReader, snps)
    return readline(vcf.v),snps
end

done(vcf::VCFReader,snps) = eof(vcf.v)

function main(args)
    filename = args[1]
    vcf = VCFReader(filename)
    snps = 0
    for (line) in vcf
        snps += 1
        if (snps %  10000 == 0)
            print(".")
        end
        if (snps % 100000 == 0)
            println(snps)
        end
    end
end

main(ARGS)