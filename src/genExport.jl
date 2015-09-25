#John Spaw
#Contains function that will convert a .dot file into a .svg image

#Converts a .dot file
function generalExport(file;filename="genImage"::String,layoutEngine="dot")
  dot = open(file,"r") do io Graph(io) end
  GraphViz.layout!(dot,engine=layoutEngine)
  open("$filename.svg","w") do f
    GraphViz.writemime(f, MIME"image/svg+xml"(),dot)
  end #do
  DEBUG && print("File saved")
end
