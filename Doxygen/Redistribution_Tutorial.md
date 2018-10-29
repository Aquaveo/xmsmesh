\tableofcontents
# Redistribution Tutorial {#Redistribution_Tutorial}

## Introduction {#Intro_Redistribution}
The purpose of this tutorial is to provide explanation on how to use the classes defined in xmsmesh to redistribute the points on a polyline according to a size function or curvature.

## Example - Redistribute Simple Polygon to Constant Size{#Example_Redistribution_Simple_Polygon}
This first example shows how to redistribute a single polygon. A picture of the example is shown below. Notice that the polygon is a simple square from (0,0) to (100,100). Also notice that the point spacing along the boundary is a constant value of 10. We will redistribute the boundary to a new constant size of 20.0.

![Simple Polygon with boundary spacing = 10.0](tutMesh_SimplePolygon_Input.png)

The basic steps to redistribute polygon points are:
1. Define the polygon as a vector of points.
2. Setup the MePolyRedistributePts class and call the Redistribute method.

\snippet xmsmesh/tutorial/TutMeshing.cpp snip_test_example_SimplePolygon_Redistribute

An image of the redistributed polygon from this example is shown below.

![Redistributed output polygon with boundary spacing = 20.0](tutRedist_SimplePolygon_Output.png)

## Example - Redistribute Simple Polygon with a Size Function{#Example_Redistribution_Size_Function}
This example show how to redistribute a polygon boundary using a size function. The size function is specified using xms::InterpBase. The InterpBase class performs spatial interpolation from points and triangles. This example uses a simple polygon with a set of 5 points and 4 triangles to define a linear size function.A picture of the example is shown below.

![Simple polygon with linear size function](tutRedist_Size_Function_Input.png)

\snippet xmsmesh/tutorial/TutMeshing.cpp snip_test_Example_Redistribute_SizeFunction

An image of the redistributed polygon from this example is shown below.

![Redistributed polygon from size function](tutRedist_Size_Function_Output.png)

## Example - Redistribute Simple Polyline using Curvature{#Example_Redistribution_Curvature}
This example shows how to redistribute a polyline using curvature. The following variables must be specified when using curvature: feature size, mean spacing, minimum curvature, and smoothing. "Feature size" is the size of the smallest feature in the polyline to be detected. "Mean spacing" is the mean spacing between the distributed points. "Minimum curvature" is the value of the curvature to be used instead of 0 in straight lines. "Smoothing" detemines if the curvatures are to be averaged by a rolling 0.25-0.5-0.25 weighted rolling average. A picture of the input polyline is shown below. We will redistribute with a feature size of 3.0, a mean spacing of 0.5, a minimum curvature of 0.0001, and no smoothing.

![Input polyline](tutRedist_Curvature_Input.png)

\snippet xmsmesh/tutorial/TutMeshing.cpp test_Example_Redistribute_Curvature

An image of the redistributed polyline from this example is shown below.

![Output polyline](tutRedist_Curvature_Output.png)

## Example - Redistribute Simple Polygon using Curvature{#Example_Redistribution_Curvature_Polygon}
This example shows how to redistribute a polygon. It is the same as the previous example except that this is a closed polygon instead of a polyline. The input is shown below.

![Input polyline](tutRedist_Curvature_Polygon_Input.png)

\snippet xmsmesh/tutorial/TutMeshing.cpp test_Example_Redistribute_Polygon_Curvature

An image of the redistributed polygon from this example is shown below.

![Output polyline](tutRedist_Curvature_Polygon_Output.png)

