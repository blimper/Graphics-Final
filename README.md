# Assignment 5: IFS Renderer

# See Project Images Here -> <http://imgur.com/a/MBm6H>

For our final project, we decided to extend/modify our SVG-drawer to render fractals from Iterated Function Systems.
The program takes in a file that includes the basic Iterated Function System as well as a text file that creates varations of said IFS.
You can also specify how many points you wish to plot. Increasing this number will make your output fractal look smoother and cleaner at the cost 
of runtime.

After parsing the IFS file, the program creates various affine transformations, each with an a attached probability.
Then, given a starting point, the program iteratively selects a random transfromation and applies it to the current point in order to plot the next.  
Without adding anything extra, this suffices for rendering a simple fractal. The steps below enhance certain aesthetic qualities of the fractals.

# Irradiance Caching:
If every plotted point is rasterized with the same intensity, the image ends up looking flat. It makes sense that pixels that are hit more often should be brighter. Irradiance caching solves this problem. Now, we initialize a buffer that stores energy bilinearly as we "plot" each point. In our resolve function, we then resample this buffer and normalize the alpha values such that we end up with varied transparencies. See image 1 on our imgur link.

# Distance Function:
In order to add some interesting color to the fractal's background, we implemented a flood-fill algorithm that colors the canvas depending on how far the pixel is away from the fractal. This algorithm loops over all the pixels and steps outward from the pixels containing fractal elements. On certain fractals, this can ceate an interesting pattern. 

# Variations:
The basic idea of the fractal flame algorithm is that after we apply the given affine functions, we run the output through another series of functions that represent variations in our fractals. These variations are predefined with varying amounts of inputs. With our given variation files, we can take a weighted sum of each of the variation functions. The files located in the "variations" folder supply the IFS renderer with information on how to weight these variations.
variation files follow this format:
```
1
0
0
0
...
```
where the nth number represents the weight of the nth variation (starting at n = 0). For a variation file to be valid, it should have k numbers, where k is the number of variations implemented (currently 16). In addition, all of the numbers should add up to 1.

# Coloring:
We assign a color to each affine function in an IFS file. To color the IFS, we calculate the new color based on the average of the color of the previous point and the color assigned to each function.

#IFS File Format
IFS files (located in the IFS folder) follow the following format:
```
n                                           (The number of affine functions)

a1 b1 c1 d1 e1 f1                           (The six coefficients of the affine function)
p1                                          (The probability of the affine function)
red1 green1 blue1                           (The RGB values assigned to the given affine function)
apost1 bpost1 cpost1 dpost1 epost1          (post transforms of the given affine function)

a2 b2 c2 d2 e2 f2
p2
red2 green2 blue2
apost2 bpost2 cpost2 dpost2 epost2

...

a b c d e f (Final Affine Transform)
```


