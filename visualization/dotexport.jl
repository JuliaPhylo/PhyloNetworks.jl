#John Spaw
#Contains function that will convert a .dot file into a .svg image

using  GraphViz
#Converts a .dot file
function dotExport(file;filename="graphimage"::String)
  #println("Opening .dot file and storing as variable")
  dot = open(file,"r") do io Graph(io) end
  #println(".dot file stored as variable")
  #println("Calling layout")
  GraphViz.layout!(dot,engine="dot")
  #println("Starting file write.")
  open("visualization/$filename.svg","w") do f
    GraphViz.writemime(f, MIME"image/svg+xml"(),dot)
  end #do
  #print("File saved")
end
