using  GraphViz
#Converts a .dot file
function dotExport(file)
  println("Opening .dot file and storing as variable")
  dot = open(file,"r") do io Graph(io) end
  println(".dot file stored as variable")
  println("Calling layout")
  GraphViz.layout!(dot,engine="dot")
  println("Starting file write.")
  open("scratchimage.svg","w") do f
    GraphViz.writemime(f, MIME"image/svg+xml"(),dot)
  end #do
  print("File saved")
end