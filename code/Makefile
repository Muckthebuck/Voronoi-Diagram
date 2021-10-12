voronoi2: voronoi2.o readData.o dcel.o output.o problem.o
	gcc -Wall -o voronoi2 voronoi2.o dcel.o readData.o output.o problem.o -g -lm

voronoi2-gui: voronoi2-gui.o readData.o dcel.o output.o gtk.o problem.o
	gcc `pkg-config --cflags gtk+-3.0` -o voronoi2-gui voronoi2-gui.o dcel.o readData.o output.o gtk.o problem.o `pkg-config --libs gtk+-3.0` -Wall -lm -g `pkg-config --static --libs --cflags igraph`

voronoi2.o: voronoi2.c watchtowerStruct.h readData.h dcel.h output.h output.h problem.h
	gcc -Wall -o voronoi2.o -c voronoi2.c -g

voronoi2-gui.o: voronoi2.c watchtowerStruct.h readData.h dcel.h output.h output.h problem.h 
	gcc -Wall -o voronoi2-gui.o -c voronoi2.c -g -DUSE_GUI

problem.o: problem.c problem.h
	gcc -Wall -o problem.o problem.c -g -c

readData.o: readData.c watchtowerStruct.c watchtowerStruct.h
	gcc -Wall -o readData.o -c readData.c -g
    
dcel.o: dcel.h dcel.c
	gcc -Wall -o dcel.o -c dcel.c -g

output.o: output.h watchtowerStruct.h dcel.h output.c
	gcc -Wall -o output.o -c output.c -g
    
gtk.o: gtk.h gtk.c
	gcc -g `pkg-config --cflags gtk+-3.0` -o gtk.o -Wall -c gtk.c `pkg-config --cflags igraph`    
