#John Spaw
#Contains function that will convert a .dot file into a .svg image

using  GraphViz
#Converts a .dot file
function generalExport(file;filename="genImage"::String,layoutEngine="dot")
  dot = open(file,"r") do io Graph(io) end
  GraphViz.layout!(dot,engine=layoutEngine)
  open("visualization/$filename.svg","w") do f
    GraphViz.writemime(f, MIME"image/svg+xml"(),dot)
  end #do
  print("File saved")
end
