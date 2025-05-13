# helpers/fileutils.jl
module FileUtils

export files_equal_chunked

function files_equal_chunked(path1::String, path2::String; chunksize=4096)::Bool
    open(path1, "r") do f1
        open(path2, "r") do f2
            while !eof(f1) && !eof(f2)
                b1 = read(f1, min(chunksize, filesize(path1)))
                b2 = read(f2, min(chunksize, filesize(path2)))
                if b1 != b2
                    return false
                end
            end
            return eof(f1) == eof(f2)
        end
    end
end

end  # module