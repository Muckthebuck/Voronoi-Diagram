/* Compile with
gcc `pkg-config --cflags gtk+-3.0` -o gtk gtk.c `pkg-config --libs gtk+-3.0` -lm -g `pkg-config --static --libs --cflags igraph`
*/

#include <gtk/gtk.h>
#include <math.h>
#include <cairo.h>
#include <igraph.h>
#include <assert.h>
#include <stdlib.h>
#include "gtk.h"

#define WIDTH   640
#define HEIGHT  480

#define ZOOM_X  100.0
#define ZOOM_Y  100.0

/* Radius of watchtower circles */
#define WTRADIUS (0.01)
/* Radius of polygon circles */
#define POLYPRADIUS (0.1)

// #define NUMVERTICES 6
// #define NUMEDGES 7
#define LAYOUTDIMENSIONS 2

/* Distance along arrow to label weight. */
#define LABEL_DISTANCE (0.5)
#define FONT_SIZE (0.2)
#define Y_OFFSET (0.05)
#define X_OFFSET (0.05)

#define VERTEX_FONT_SIZE (0.2)

gfloat f (gfloat x)
{
    return 0.03 * pow (x, 3);
}

igraph_t graph;
igraph_vector_t edges;
// Edges in shortest path 0 excluded, 1 included.
igraph_vector_t shortest_path;
igraph_matrix_t res;

gdouble minx = -1;
gdouble miny = -1;
gdouble maxx = 1;
gdouble maxy = 1;

int selectedWT = -1;
int selectedPolyPoint = -1;
int selectedPolyEdge = -1;

double pointer_x = 0;
double pointer_y = 0;

int gtkEdgeCount = 0;
int gtkNodeCount = 0;

void **dijkstraResults = NULL;

// GtkWidget **weightLabels = NULL;
char **weightStrings = NULL;
char **gtkVertexLabels = NULL;

void setupGTKBounds(double newxmin, double newxmax, double newymin, double newymax){
    minx = newxmin;
    miny = newymin;
    maxx = newxmax;
    maxy = newymax;
}

/* Watchtower point data and fetch functions. */
int pointsSet = 0;
int gtkpointCount = 0;
void *gtkdisplayData = NULL;
char *(*gtkgetDataString)(void *, int) = NULL;
void (*gtkfreeDataString)(char *) = NULL;
double (*gtkgetPointx)(int, void *) = NULL;
double (*gtkgetPointy)(int, void *) = NULL;

/* Polygon points. */
int polygonsPointsSet = 0;
int gtkdcelVerticesCount;
void *gtkdcel;
double (*gtkgetDCELVertexX)(void *, int);
double (*gtkgetDCELVertexY)(void *, int);

/* Polygon lines. */
int polygonLinesSet = 0;
int gtkedgeCount;
void *gtkedgedcel;
int (*gtkgetDCELEdgeVertexStart)(void *, int);
int (*gtkgetDCELEdgeVertexPairStart)(void *, int);
int (*gtkgetDCELEdgeVertexEnd)(void *, int);
int (*gtkgetDCELEdgeVertexPairEnd)(void *, int);
int (*gtkDCELhasEdge)(void *, int);
int (*gtkDCELhasEdgePair)(void *, int);

/* GTK Action & context */
int actionSet = 0;
void (*gtkdoAction)(void *, double, double);
void *gtkactionContext;

void setGTKWatchTower(int vertex){
    selectedWT = vertex;
}

void setGTKVertex(int vertex){
    selectedPolyPoint = vertex;
}

void setGTKEdge(int edge){
    selectedPolyEdge = edge;
}

void setupGTKPoints(int pointCount, void *displayData, char *(*getDataString)(void *, int), 
    void (*freeDataString)(char *), double (*getPointx)(int, void *), 
    double (*getPointy)(int, void *)){
    pointsSet = 1;
    gtkpointCount = pointCount;
    gtkdisplayData = displayData;
    gtkgetDataString = getDataString;
    gtkfreeDataString = freeDataString;
    gtkgetPointx = getPointx;
    gtkgetPointy = getPointy;
}

void setupGTKPolyPoints(int dcelVerticesCount, void *dcel, 
    double (*getDCELVertexX)(void *, int), double (*getDCELVertexY)(void *, int)){
    polygonsPointsSet = 1;
    gtkdcelVerticesCount = dcelVerticesCount;
    gtkdcel = dcel;
    gtkgetDCELVertexX = getDCELVertexX;
    gtkgetDCELVertexY = getDCELVertexY;
}

void setupGTKPolyLines(int edgeCount, void *dcel, 
    int (*getDCELEdgeVertexStart)(void *, int),
    int (*getDCELEdgeVertexPairStart)(void *, int),
    int (*getDCELEdgeVertexEnd)(void *, int),
    int (*getDCELEdgeVertexPairEnd)(void *, int), int (*DCELhasEdge)(void *, int),
    int (*DCELhasEdgePair)(void *, int)){
    polygonLinesSet = 1;
    gtkedgeCount = edgeCount;
    gtkedgedcel = dcel;
    gtkgetDCELEdgeVertexStart = getDCELEdgeVertexStart;
    gtkgetDCELEdgeVertexPairStart = getDCELEdgeVertexPairStart;
    gtkgetDCELEdgeVertexEnd = getDCELEdgeVertexEnd;
    gtkgetDCELEdgeVertexPairEnd = getDCELEdgeVertexPairEnd;
    gtkDCELhasEdge = DCELhasEdge;
    gtkDCELhasEdgePair = DCELhasEdgePair;
}

void setupGTKAction(void (*doAction)(void *, double, double), void *actionContext){
    actionSet = 1;
    gtkdoAction = doAction;
    gtkactionContext = actionContext;
}

void update_destination(){
    /* TODO: Implement hover mode handling. */
    return;
}

static gboolean
on_draw (GtkWidget *widget, cairo_t *cr, gpointer user_data)
{
    GdkRectangle da;            /* GtkDrawingArea size */
    gdouble dx = 5.0, dy = 5.0; /* Pixels between each point */
    gdouble i, clip_x1 = 0.0, clip_y1 = 0.0, clip_x2 = 0.0, clip_y2 = 0.0;
    GdkWindow *window = gtk_widget_get_window(widget);

    /* Determine GtkDrawingArea dimensions */
    gdk_window_get_geometry (window,
            &da.x,
            &da.y,
            &da.width,
            &da.height);
    
    /***************************************/

    /* Draw on a black background */
    cairo_set_source_rgb (cr, 0.0, 0.0, 0.0);
    cairo_paint (cr);

    /* Change the transformation matrix */
    cairo_translate (cr, da.width / 2, da.height / 2);
    cairo_scale (cr, ZOOM_X, ZOOM_Y);

    /* Determine the data points to calculate (ie. those in the clipping zone) */
    cairo_device_to_user_distance (cr, &dx, &dy);
    cairo_clip_extents (cr, &clip_x1, &clip_y1, &clip_x2, &clip_y2);
    cairo_set_line_width (cr, dx);

    /* Draw watchtower points */
    double x, y;
    int j, k;
    double rangeX = maxx - minx;
    double rangeY = maxy - miny;
    double avgX = (maxx + minx)/2;
    double avgY = (maxy + miny)/2;
    if(pointsSet){
        for(j = 0; j < gtkpointCount; j++){
            x = ((gtkgetPointx(j, gtkdisplayData) - avgX) / rangeX) * 5;
            y = (-(gtkgetPointy(j, gtkdisplayData) - avgY) / rangeY) * 5;
            cairo_move_to (cr, x, y);
            cairo_new_sub_path(cr);
            // Draw selected source in yellow.
            if (j == selectedWT){
                cairo_set_source_rgba (cr, 1, 0.6, 0.2, 0.6);
            } else {
                cairo_set_source_rgba (cr, 1, 0.2, 0.2, 0.6);
            }

            cairo_arc(cr, x, y, WTRADIUS, 0, 2 * M_PI);
            cairo_stroke(cr);
        }
    }
    
    /* Draw points for polygon */
    if(polygonsPointsSet){
        for(j = 0; j < gtkdcelVerticesCount; j++){
            x = ((gtkgetDCELVertexX(gtkdcel, j) - avgX) / rangeX) * 5;
            y = (-(gtkgetDCELVertexY(gtkdcel, j) - avgY) / rangeY) * 5;
            cairo_move_to (cr, x, y);
            cairo_new_sub_path(cr);
            // Draw selected source in yellow.
            if (j == selectedPolyPoint){
                cairo_set_source_rgba (cr, 1, 0.6, 0.2, 0.6);
            } else {
                cairo_set_source_rgba (cr, 0.2, 0.6, 0.6, 0.6);
            }

            cairo_arc(cr, x, y, POLYPRADIUS, 0, 2 * M_PI);
            cairo_stroke(cr);
        }
    }
    
    /* Draw lines for polygon */
    int startVertex;
    int endVertex;
    gdouble xStart;
    gdouble yStart;
    gdouble xEnd;
    gdouble yEnd;
    int pair;
    int pairStartVertex;
    int pairEndVertex;
    gdouble pairxStart;
    gdouble pairyStart;
    gdouble pairxEnd;
    gdouble pairyEnd;
    
    /* Immediate drawing context */
    //gdouble x_from, x_to, y_from, y_to;
    void drawLine(gdouble x_from, gdouble y_from, gdouble x_to, gdouble y_to, int edge,
        int sVertex, int eVertex){
        // Draw arrowhead code adapted from 
        // http://kapo-cpp.blogspot.com/2008/10/drawing-arrows-with-cairo.html
        cairo_move_to (cr, (x_from), (y_from));
        if(sVertex == selectedPolyPoint || eVertex == selectedPolyPoint || 
           edge == selectedPolyEdge){
            cairo_set_source_rgba (cr, 0.6, 0.6, 0.6, 0.25);
        } else {
            cairo_set_source_rgba (cr, 0.2, 0.6, 0.6, 0.25);
        }
        double angle = atan2 (y_to - y_from, x_to - x_from) + M_PI;
        cairo_line_to (cr, (x_to + (POLYPRADIUS / 2.0) * cos(angle)), 
                           (y_to + (POLYPRADIUS / 2.0) * sin(angle)));
        cairo_stroke(cr);

        double arrow_degrees = 0.5;
        double x1 = x_to + 0.05 * cos(angle) + 0.2 * cos(angle - arrow_degrees);
        double y1 = y_to + 0.05 * sin(angle) + 0.2 * sin(angle - arrow_degrees);
        double x2 = x_to + 0.05 * cos(angle) + 0.2 * cos(angle + arrow_degrees);
        double y2 = y_to + 0.05 * sin(angle) + 0.2 * sin(angle + arrow_degrees);
        cairo_move_to (cr, x_to + 0.05 * cos(angle), y_to + 0.05 * sin(angle));
        cairo_line_to (cr, x1, y1);
        cairo_stroke(cr);
        
        /* Only draw one side of line head. */
        //cairo_move_to (cr, x_to + 0.05 * cos(angle), y_to + 0.05 * sin(angle));
        //cairo_line_to (cr, x2, y2);
        //cairo_stroke(cr);
    }
    if(polygonLinesSet){
        if(! polygonsPointsSet){
            fprintf(stderr, "Polygon points not set, can't construct edges\n");
            assert(polygonsPointsSet);
        }
        for(j = 0; j < gtkedgeCount; j++){
            if(! gtkDCELhasEdge(gtkedgedcel, j)){
                continue;
            }
            startVertex = gtkgetDCELEdgeVertexStart(gtkedgedcel, j);
            endVertex = gtkgetDCELEdgeVertexEnd(gtkedgedcel, j);
            xStart = ((gtkgetDCELVertexX(gtkdcel, startVertex) - avgX) / rangeX) * 5;
            yStart = (-(gtkgetDCELVertexY(gtkdcel, startVertex) - avgY) / rangeY) * 5;
            xEnd = ((gtkgetDCELVertexX(gtkdcel, endVertex) - avgX) / rangeX) * 5;
            yEnd = (-(gtkgetDCELVertexY(gtkdcel, endVertex) - avgY) / rangeY) * 5;
            pair = gtkDCELhasEdgePair(gtkedgedcel, j);
            if(pair){
                pairStartVertex = gtkgetDCELEdgeVertexPairStart(gtkedgedcel, j);
                pairEndVertex = gtkgetDCELEdgeVertexPairEnd(gtkedgedcel, j);
                pairxStart = 
                    ((gtkgetDCELVertexX(gtkdcel, pairStartVertex) - avgX) / rangeX) * 5;
                pairyStart = 
                    (-(gtkgetDCELVertexY(gtkdcel, pairStartVertex) - avgY) / rangeY) * 5;
                pairxEnd = 
                    ((gtkgetDCELVertexX(gtkdcel, pairEndVertex) - avgX) / rangeX) * 5;;
                pairyEnd = 
                    (-(gtkgetDCELVertexY(gtkdcel, pairEndVertex) - avgY) / rangeY) * 5;
            }
            
            /* Draw lines */
            
            drawLine(xStart, yStart, xEnd, yEnd, j, startVertex, endVertex);
            if(pair){
                drawLine(pairxStart, pairyStart, pairxEnd, pairyEnd, j, 
                         pairStartVertex, pairEndVertex);
            }
        }
    }
    
    /* Draw cursor circle */
    cairo_move_to (cr, pointer_x, pointer_y);
    cairo_new_sub_path(cr);
    cairo_set_source_rgba (cr, 0.2, 0.2, 0.2, 0.6);
    cairo_arc(cr, pointer_x, pointer_y, 0.1, 0, M_PI / 2);
    cairo_stroke(cr);
    cairo_arc(cr, pointer_x, pointer_y, 0.1, M_PI, 3 * M_PI / 2);
    cairo_stroke(cr);

    return FALSE;
}

static gboolean mouse_moved(GtkWidget *widget, GdkEvent *event, gpointer user_data) {

    if (event->type==GDK_MOTION_NOTIFY) {
        GdkEventMotion* e=(GdkEventMotion*)event;
        // printf("Coordinates: (%f,%f)\n", e->x, e->y);
        GdkRectangle da;            /* GtkDrawingArea size */
        gdouble dx = 5.0, dy = 5.0; /* Pixels between each point */
        gdouble i, clip_x1 = 0.0, clip_y1 = 0.0, clip_x2 = 0.0, clip_y2 = 0.0;
        GdkWindow *window = gtk_widget_get_window(widget);

        /* Determine GtkDrawingArea dimensions */
        gdk_window_get_geometry (window,
                &da.x,
                &da.y,
                &da.width,
                &da.height);
        // printf("%f %f %d %d\n", (e->x - da.width / 2) / dx, (e->y - da.height / 2) / dy, da.width, da.height);
        // Scale to Cairo zoom level
        pointer_x = (e->x - da.width / 2) / ZOOM_X;
        pointer_y = (e->y - da.height / 2) / (ZOOM_Y);
        update_destination();
        
        gtk_widget_queue_draw(widget);
    }
    return FALSE;
}

static gboolean button_pressed(GtkWidget *widget, GdkEvent *event, gpointer user_data){
    if(actionSet){
        double xrange = maxx - minx;
        double xavg = (minx + maxx)/2;
        double yrange = maxy - miny;
        double yavg = (miny + maxy)/2;
        double transformedX = pointer_x / 5 * xrange + xavg;
        double transformedY = -pointer_y / 5 * yrange + yavg;
        gtkdoAction(gtkactionContext, transformedX, transformedY);
    }
    
    return FALSE;
}

GtkWidget *window = NULL;
GtkWidget *da = NULL;

void setupGTK(int *argc, char ***argv, char *title){
    gtk_init (argc, argv);
    window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_default_size(GTK_WINDOW (window), WIDTH, HEIGHT);
    gtk_window_set_title(GTK_WINDOW (window), title);

    g_signal_connect(G_OBJECT (window), "destroy", gtk_main_quit, NULL);

    da = gtk_drawing_area_new();

    gtk_container_add(GTK_CONTAINER(window), da);
}

void runGTKInteraction(){
    g_signal_connect (G_OBJECT (da),
            "draw",
            G_CALLBACK (on_draw),
            NULL);
    
    g_signal_connect (G_OBJECT (da),
            "motion-notify-event",
            G_CALLBACK (mouse_moved),
            NULL);

    g_signal_connect (G_OBJECT (window),
            "button-release-event",
            G_CALLBACK (button_pressed),
            NULL);
    
    gtk_widget_set_events(da, GDK_POINTER_MOTION_MASK);
    gtk_widget_add_events(window, GDK_BUTTON_RELEASE_MASK);
    gtk_widget_show_all(window);
    gtk_main();
}

void freeGTK(){

}
