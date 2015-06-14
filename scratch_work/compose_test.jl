using Compose

a = compose(context(), fill("tomato"),
	(context(0.0, 0.0, 0.5, 0.5), circle()),
	(context(0.0, 1.0, 0.1, 0.1), circle()),
	(context(1.0, 0.0, 0.3, 0.3), circle()))

#a = compose(context(), fill("tomato"),
#        (context(0.0, 0.2, 0.5, 0.5), circle()),
#        (context(0.5, 0.5, 0.5, 0.5), circle()),
#        (context(1.0, 0.0, 0.3, 0.6), rectangle()))

img = SVG("test_a.svg", 10inch, 10inch)
draw(img,a)

b = compose(context(), fill("bisque"),
	(context(3.0, 0.0, 0.5, 0.5), rectangle()))
	
draw(img,a+b)