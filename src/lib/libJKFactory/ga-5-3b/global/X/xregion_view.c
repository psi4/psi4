#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include "xregion.h"

void create_main_window()
{
  Display *disp;
  int screen;

  disp = XtDisplay(overview_shell);
  screen = DefaultScreen(disp);
 
  /* Create top level shell widget */
  top_level_widget = XtVaAppCreateShell("xregion","XRegion",
                     applicationShellWidgetClass,disp,
                     NULL);

  /* Create form widget to hold everything else */

  box_widget = XtVaCreateManagedWidget("box", formWidgetClass,
                                     top_level_widget,
                                     XtNbackground, DEFAULT_BG,
                                     XtNforeground, DEFAULT_FG,
                                     NULL);

  /* Create the label to hold the title */

  (void) strcpy(title, "Array Access Display");
  title_widget = XtVaCreateManagedWidget("title", labelWidgetClass,
                                       box_widget,
                                       XtNbackground, DEFAULT_BG,
                                       XtNforeground, DEFAULT_FG,
                                       XtNx, 10,
                                       XtNy, 5,
                                       XtNwidth, 300,
                                       XtNlabel, title,
                                       XtNborderWidth, 0,
                                       NULL);

  coord_widget = XtVaCreateManagedWidget("coords", labelWidgetClass,
                                       box_widget,
                                       XtNbackground, DEFAULT_BG,
                                       XtNforeground, DEFAULT_FG,
                                       XtNy, 5,
                                       XtNhorizDistance, 10,
                                       XtNfromHoriz, title_widget,
                                       XtNjustify, XtJustifyLeft,
                                       XtNlabel, "Coordinates x, y:                 ",
                                       XtNborderWidth, 0,
                                       NULL);

  /* Create the Quit command button */

  quit_button = XtVaCreateManagedWidget("quit", commandWidgetClass,
                                       box_widget,
                                       XtNbackground, DEFAULT_BG,
                                       XtNforeground, DEFAULT_FG,
                                       XtNx, 10,
                                       XtNvertDistance, 15,
                                       XtNfromVert, title_widget,
                                       XtNlabel, "Quit",
                                       XtNshapeStyle, XmuShapeOval,
                                       NULL);
  XtAddCallback(quit_button, XtNcallback, Quit, NULL);

  /* Create the Start/Stop command button */

  start_stop_button = XtVaCreateManagedWidget("start/stop", commandWidgetClass,
                                       box_widget,
                                       XtNbackground, DEFAULT_BG,
                                       XtNforeground, DEFAULT_FG,
                                       XtNvertDistance, 15,
                                       XtNfromVert, title_widget,
                                       XtNhorizDistance, 10,
                                       XtNfromHoriz, quit_button,
                                       XtNlabel, "Start",
                                       XtNshapeStyle, XmuShapeOval,
                                       NULL);
  XtAddCallback(start_stop_button, XtNcallback, StartStop, NULL);

  /* Create the scroll bar for the interval */

  interval_label = XtVaCreateManagedWidget("interval_label", labelWidgetClass,
                                       box_widget,
                                       XtNbackground, DEFAULT_BG,
                                       XtNforeground, DEFAULT_FG,
                                       XtNvertDistance, 5,
                                       XtNfromVert, title_widget,
                                       XtNhorizDistance, 20,
                                       XtNfromHoriz, start_stop_button,
                                       XtNlabel, "Time Interval",
                                       XtNborderWidth, 0,
                                       NULL);

  scroll_widget = XtVaCreateManagedWidget("scroll", scrollbarWidgetClass,
                                        box_widget,
                                       XtNbackground, DEFAULT_BG,
                                       XtNforeground, DEFAULT_FG,
                                       XtNvertDistance, 5,
                                       XtNfromVert, interval_label,
                                       XtNhorizDistance, 20,
                                       XtNfromHoriz, start_stop_button,
                                       XtNorientation, XtorientHorizontal,
                                       XtNlength, 100,
                                       XtNthickness, 15,
                                       NULL);
  XtAddCallback(scroll_widget, XtNscrollProc, ScrollProc, NULL);
  XtAddCallback(scroll_widget, XtNjumpProc, JumpProc, NULL);

  /* Create the label widget which displays the interval value
     associated with the scrollbar. */

  (void) sprintf(interval_string, "%4d ms", interval);
  interval_widget = XtVaCreateManagedWidget("interval", labelWidgetClass,
                                           box_widget,
                                       XtNbackground, DEFAULT_BG,
                                       XtNforeground, DEFAULT_FG,
                                       XtNvertDistance, 5,
                                       XtNfromVert, interval_label,
                                       XtNhorizDistance, 5,
                                       XtNfromHoriz, scroll_widget,
                                       XtNjustify, XtJustifyRight,
                                       XtNlabel, interval_string,
                                       XtNborderWidth, 0,
                                       NULL);

  /* Create the scroll bar for the slowdown */

  slowdown_label = XtVaCreateManagedWidget("slowdown_label", labelWidgetClass,
                                       box_widget,
                                       XtNbackground, DEFAULT_BG,
                                       XtNforeground, DEFAULT_FG,
                                       XtNvertDistance, 5,
                                       XtNfromVert, title_widget,
                                       XtNhorizDistance, 25,
                                       XtNfromHoriz, interval_widget,
                                       XtNlabel, "Slowdown Factor",
                                       XtNborderWidth, 0,
                                       NULL);

  scroll_widget2 = XtVaCreateManagedWidget("scroll2", scrollbarWidgetClass,
                                        box_widget,
                                       XtNbackground, DEFAULT_BG,
                                       XtNforeground, DEFAULT_FG,
                                       XtNvertDistance, 5,
                                       XtNfromVert, slowdown_label,
                                       XtNhorizDistance, 25,
                                       XtNfromHoriz, interval_widget,
                                       XtNorientation, XtorientHorizontal,
                                       XtNlength, 100,
                                       XtNthickness, 15,
                                       NULL);
  XtAddCallback(scroll_widget2, XtNscrollProc, ScrollProc2, NULL);
  XtAddCallback(scroll_widget2, XtNjumpProc, JumpProc2, NULL);

  /* Create the label widget which displays the slowdown value
     associated with the scrollbar 2. */

  (void)  sprintf(slowdown_string, "%5d times", (long)slowdown);
  slowdown_widget = XtVaCreateManagedWidget("slowdown", labelWidgetClass,
                                           box_widget,
                                       XtNbackground, DEFAULT_BG,
                                       XtNforeground, DEFAULT_FG,
                                       XtNvertDistance, 5,
                                       XtNfromVert, slowdown_label,
                                       XtNhorizDistance, 5,
                                       XtNfromHoriz, scroll_widget2,
                                       XtNjustify, XtJustifyRight,
                                       XtNlabel, slowdown_string,
                                       XtNborderWidth, 0,
                                       NULL);

  /* Now add the actual canvas ... */

  canvas_widget = XtVaCreateManagedWidget("canvas", formWidgetClass,
                                        box_widget,
                                        XtNheight, pict_height * scale,
                                        XtNwidth, pict_width * scale, 
                                        XtNvertDistance, 20,
                                        XtNfromVert, quit_button,
                                        XtNbackground, CANVAS_COLOR,
                                        XtNborderWidth, 0,
                                        NULL);
  /* Add callback for exposure */

  XtAddEventHandler(canvas_widget,ExposureMask,False,Exposed,NULL);
  XtAddEventHandler(canvas_widget,PointerMotionMask,False,running_coords,NULL);

  /* Now add the color scale ... */

  map_widget = XtVaCreateManagedWidget("colorMap", compositeWidgetClass,
                                        box_widget,
                                       XtNbackground, DEFAULT_BG,
                                       XtNforeground, DEFAULT_FG,
                                       XtNheight, 350,
                                       XtNwidth, 80, 
                                       XtNvertDistance, 20,
                                       XtNfromVert, quit_button,
                                       XtNhorizDistance, 20,
                                       XtNfromHoriz, canvas_widget,
                                       XtNborderWidth, 0,
                                       NULL);

}

/* JJU: void running_coords(Widget widget, caddr_t data, XEvent *event) */
void running_coords(Widget widget, XtPointer data, XEvent *event, Boolean *bln)
{
  char loc_str[40];
  int x, y;

  x = event->xmotion.x;
  y = event->xmotion.y;

/*   sprintf(loc_str, "Coordinates x, y: %d, %d", (int) (x / scale) + left_edge + 1, 
         (int) (y / scale) + top_edge + 1);
*/
  sprintf(loc_str, "Coordinates x, y: %d, %d", 
		(int) (x / scale) + left_edge, 
         	(int) (y / scale) + top_edge);
  XtVaSetValues(coord_widget, XtNlabel, loc_str, NULL);
}

void setup_drawing()
{
  int i, x, y;
  XGCValues gcv;

  /* Set up the drawing environment */

  display = XtDisplay(canvas_widget);
  window = XtWindow(canvas_widget);
  window_map = XtWindow(map_widget);
  screen = DefaultScreen(display);
  visual = DefaultVisual(display, screen);
  depth = DisplayPlanes(display, screen);
  (void) printf("depth = %d\n",depth);

  gc = XCreateGC(display, window, 0, (XGCValues *) NULL);
  XSetForeground(display, gc, GRID_COLOR);

  gcv.font = XLoadFont(display, "8x13");
  if(!gcv.font)
  {
    printf("error font not loaded\n");
  }

  gc_map = XCreateGC(display, window_map, GCFont, &gcv);

  Setcmap();

  /* Make image to match the size of our canvas */
  x = pict_width * scale;
  y = pict_height * scale;
  pict  = (u_char *) malloc((x + pict_width) * y);
  image = XCreateImage(display, visual, depth, ZPixmap, 0,
                       pict, x, y, 8, x);

  /* Make the byte array which will hold the access data */

  if (!(grid = (u_char *) malloc((unsigned) (pict_width * pict_height))))
  {
    Error("failed to allocate grid", -1);
  }
  bzero((char *) grid, pict_width * pict_height);

  /* Make the byte array which will hold the access flag */

  if (!(flag = (u_char *) malloc((unsigned) (pict_width * pict_height))))
  {
    Error("failed to allocate flag", -1);
  }
  bzero((char *) flag, pict_width * pict_height);

  /* Make the array which will hold the integral */

  if (!(integr = (double *) malloc(sizeof(double) * (pict_width * pict_height))))
  {
    Error("failed to allocate integr", -1);
  }

  /* Make the array which will hold the last access time */

  if (!(ltime = (double *) malloc(sizeof(double) * (pict_width * pict_height))))
  {
    Error("failed to allocate ltime", -1);
  }

  for(i = 0; i < pict_width * pict_height; i++, *ltime = 0.0, *integr = 0.0);

  /* clear the array display */
  UpdatePixRegion(0, pict_height - 1, 0, pict_width - 1, 0, 0.0);
  DisplayPixRegion(0, pict_height - 1, 0, pict_width - 1);
}

/**/
/* JJU: void Quit(Widget widget, caddr_t data, XEvent *event) */
void Quit(Widget widget, XtPointer data, XtPointer event)
{
  exit(0);
}


/**/
/* JJU: void StartStop(Widget widget, caddr_t data, XEvent *event) */
void StartStop(Widget widget, XtPointer data, XtPointer event)
{
  /* Toggle propagation of display */

  if (working)
  {
    XtRemoveTimeOut(timer);
    working = False;
    XtSetArg(arg[0], XtNlabel, "Start");   /* Reset button label */
    XtSetValues(start_stop_button,arg,1);
    XFlush(display);
  }
  else
  {
    XtSetArg(arg[0], XtNlabel, "Stop");   /* Reset button label */
    XtSetValues(start_stop_button,arg,1);
    timer = XtAppAddTimeOut(xregion_app, interval, TimeOutCallback, NULL);
    working = True;
    XFlush(display);
  }
}

