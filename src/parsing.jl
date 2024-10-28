using XLSX

function parseMailXLSX(path::String)
    file_xlsx = XLSX.readxlsx(path)

    # Getting the names of all sheets in the file
    names = XLSX.sheetnames(file_xlsx)

    for e in names
        println(e)
    end

    stringParameters = split(names[1], '_')
    parameters = [parse(Int64,e) for e in stringParameters]
    cols = parameters[1]
    rows = parameters[2]

    parameters

    sh = file_xlsx[names[1]]

    mat::Matrix{Int64} = zeros(rows,cols)
    for row in 2:rows+1
        for col in 2:cols+1
            mat[row-1,col-1] = sh[row, col]
        end
    end

    return mat
end

# function parseAll(directoryPath::String)
#     all_files = readdir(directoryPath)
#     for e in all_files
#         full_path = directoryPath * e
#         parseMailXLSX(full_path)
#     end
#     #TODO: store all matrices with a dictionnary
# end