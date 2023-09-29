using Plots
using DelimitedFiles
using LaTeXStrings
using ReadableRegex
using SparseArrays

plotlyjs()

function extractIntegersFromString(str)
     # Remove any non-digit characters from the string
     cleaned_str = filter(x -> isdigit(x) || x == ',', str)

     # Split the cleaned string by commas
     values = split(cleaned_str, ',')

     # Convert the values to integers
     integers = parse.(Int, values)

     return integers
end



filelist = readdir("/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs/system/")
plot()
for file in filelist
     str = "/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs/system/" * file
     if (occursin("matrix", str))
          local ff = readdlm(str)
          M=zeros(size(ff)[1],size(ff)[1])
          for i in 1:size(ff)[1]     
               pos = ff[i, 1]
               val = ff[i, 2]
               pos = String(pos)
               integers = extractIntegersFromString(pos)
               M[integers[1]+1,integers[2]+1]= val
               
          end
          sparse_matrix = sparse(M)
          rows, cols, values = findnz(sparse_matrix)
          scatter(cols, rows, zcolor=values, markersize=2 ,legend=false,colorbar = true,mode="markers")
          #heatmap(cols,rows,values)
          #heatmap(M)
          xlabel!("Columns")
          ylabel!("Rows")
          title!("Sparsity Pattern")
     end
end

current()




