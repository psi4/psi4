#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include "xregion.h"

XtAppContext create_overview(int argc, char **argv)
{
  XtAppContext app;

  overview_shell = XtVaAppInitialize(&app,"XRegion", NULL, 0, &argc, argv,
                                       NULL, XtNtitle, "select region to view",
                                       NULL);
  
  set_config();

  overview_widget = XtVaCreateManagedWidget("selectCanvas", formWidgetClass,
                                        overview_shell,
                                        XtNbackground, DEFAULT_BG,
                                        XtNforeground, DEFAULT_FG,
                                        NULL);
  
  overview_title = XtVaCreateManagedWidget("selectTitle",labelWidgetClass,
                                        overview_widget, 
                                        XtNlabel, "x, y:",
                                        XtNbackground, DEFAULT_BG,
                                        XtNforeground, DEFAULT_FG,
                                        XtNx, 10,
                                        XtNy, 10,
                                        XtNwidth, 350,
                                        XtNvertDistance, 10,
                                        XtNhorizDistance, 10,
                                        XtNborderWidth, 0,
                                        NULL);

  /* Create the Start/Stop command button */

  view_button = XtVaCreateManagedWidget("viewbutton", commandWidgetClass,
                                       overview_widget,
                                        XtNbackground, DEFAULT_BG,
                                        XtNforeground, DEFAULT_FG,
                                       XtNhorizDistance, 10,
                                       XtNvertDistance, 10,
                                       XtNfromHoriz, overview_title,
                                       XtNlabel, "View Selection",
                                       XtNsensitive, False,
                                       XtNshapeStyle, XmuShapeOval,
                                       NULL);
  XtAddCallback(view_button, XtNcallback, start_view, NULL);

  select_widget = XtVaCreateManagedWidget("areaSelect", formWidgetClass,
                                        overview_widget, 
                                        XtNvertDistance, 50,
                                        XtNhorizDistance, 5,
                                        XtNheight, overview_height,
                                        XtNwidth, overview_width,
                                        XtNborderWidth, 0,
                                        NULL);

  XtAddEventHandler(select_widget, PointerMotionMask, False, running_overview, NULL);
  XtAddEventHandler(select_widget, ButtonPressMask, False, draw_select_box, NULL);
  XtAddEventHandler(select_widget, ButtonReleaseMask, False, draw_select_box, NULL);
  XtAddEventHandler(select_widget, Button1MotionMask, False, draw_select_box, NULL);
  XtAddEventHandler(select_widget, ExposureMask, False, draw_select_box, NULL);

  return(app);
}

/* JJU: void start_view(Widget widget, caddr_t data, XEvent *event); */
void start_view(Widget widget, XtPointer data, XtPointer event)
{
  create_main_window();

  /* Realize everything */
  XtRealizeWidget(top_level_widget);

  setup_drawing();

  XtRemoveEventHandler(select_widget, ButtonPressMask, False, 
                               draw_select_box, NULL);
  XtRemoveEventHandler(select_widget, ButtonReleaseMask, False, 
                               draw_select_box, NULL);
  XtRemoveEventHandler(select_widget, Button1MotionMask, False, 
                               draw_select_box, NULL);
}

/* JJU: void running_overview(Widget widget, caddr_t data, XEvent *event) */
void running_overview(Widget widget, XtPointer data, XEvent *event,
                      Boolean *bln)
{
  char loc_str[40];
  int x, y;
  int scale_x, scale_y;
 
  x = event->xmotion.x;
  y = event->xmotion.y;

  if (x < 0) x = 0;
  if (x > overview_width - 1) x = overview_width - 1;
  if (y < 0) y = 0;
  if (y > overview_height - 1) y = overview_height - 1;
/*
  scale_x = (int)((x / overview_scale) + 1.5); 
  scale_y = (int)((y / overview_scale) + 1.5); 
*/
  scale_x = (int)((x / overview_scale) + .5); 
  scale_y = (int)((y / overview_scale) + .5); 

  if (scale_x > grid_x) scale_x = grid_x;
  if (scale_y > grid_y) scale_y = grid_y;

  sprintf(loc_str, "x, y: %d, %d", scale_x, scale_y); 

  XtVaSetValues(overview_title, XtNlabel, loc_str, NULL);
}

/* JJU: void draw_select_box(Widget widget, caddr_t data, XEvent *event) */
void draw_select_box(Widget widget, XtPointer data, XEvent *event,
                     Boolean *bln)
{
  static int x1, x2, y1, y2;
  static int first = True;
  static GC gc_rband, gc_select;
  static Display *disp;
  static Pixmap pixmap;
  /* JJU: static Window *win; */
  static Window win;
  static int screen;
  static int button_state = 0;
  static int another_state = False;
  char err_str[80];
  char loc_str[40];
  int i, ul_x, ul_y, lr_x, lr_y; 

  if (first)
  {
    first = False;
    x1 = x2 = y1 = y2 = -1;

    disp = XtDisplay(select_widget);
    win = XtWindow(select_widget);
    screen = DefaultScreen(disp);

    gc_rband = XCreateGC(disp, win, NULL, NULL);
    gc_select = XCreateGC(disp, win, NULL, NULL);

    pixmap = XCreatePixmap(disp, RootWindow(disp, screen), overview_width, 
                           overview_height, DefaultDepth(disp, screen));
    XSetForeground(disp, gc_select, CANVAS_COLOR);
    XFillRectangle(disp, pixmap, gc_select, 0, 0, overview_width, 
                   overview_height);
    XSetForeground(disp, gc_rband, RBAND_COLOR);
    XSetForeground(disp, gc_select, SELECT_COLOR);
    XSetBackground(disp, gc_rband, CANVAS_COLOR);
    XSetFunction(disp,gc_rband, GXxor);
/*
    XSetLineAttributes(disp, gc_rband, 3, LineSolid, CapRound, JoinRound); 
    XSetLineAttributes(disp, gc_select, 3, LineSolid, CapRound, JoinRound); 
*/
  }

  switch (event->type)
  {
    case Expose:
      if (event->xexpose.count == 0)
      {
        XSetForeground(disp, gc_select, GRID_COLOR);
        for(i = 0; i < rows; i++)
        {
          int _x1=0,_x2=overview_width - 1;
          int _y1,_y2;
          _y1=_y2=overlay_row[i] * overview_scale;
          XDrawLine(disp, pixmap, gc_select, _x1, _y1, _x2, _y2);
        }
        for(i = 0; i < cols; i++)
        {
          int _y1=0,_y2=overview_height - 1;
          int _x1,_x2;
          _x1 =_x2 =overlay_col[i] * overview_scale;
          XDrawLine(disp, pixmap, gc_select, _x1, _y1, _x2, _y2);
        }
        XSetForeground(disp, gc_select, SELECT_COLOR);
        XCopyArea(disp, pixmap, win, gc_select, 0, 0, overview_width, 
                  overview_height, 0, 0); 
      }
      break;

    case ButtonPress:
      if (event->xbutton.button == 1) 
      {
        x1 = x2 = event->xbutton.x;
        y1 = y2 = event->xbutton.y;
        button_state = 1;
        another_state = True;
      }
      break;

    case ButtonRelease:
      if (event->xbutton.button == 1 && button_state == 1) 
      {
        /* erase old rubberband box */
        XDrawRectangle(disp, win, gc_rband, x1, y1, x2 - x1, y2 - y1); 

        /* get latest corner */
        x2 = event->xbutton.x;
        y2 = event->xbutton.y;

        if (x2 < 0) x2 = 0;
        if (x2 > overview_width - 1) x2 = overview_width - 1;
        if (y2 < 0) y2 = 0;
        if (y2 > overview_height - 1) y2 = overview_height - 1;

        if ((abs(x2 - x1) > 4) || (abs(y2 - y1) > 4))
        {
          XtVaSetValues(overview_title, XtNlabel, loc_str, NULL);

          /* draw final selection box */
          XDrawRectangle(disp, win, gc_select, x1, y1, x2 - x1, y2 - y1); 
          XDrawRectangle(disp, pixmap, gc_select, x1, y1, x2 - x1, y2 - y1); 
          XFillRectangle(disp, pixmap, gc_select, x1, y1, x2 - x1, y2 - y1); 
          XSetForeground(disp, gc_select, GRID_COLOR);
        for(i = 0; i < rows; i++)
        {
          int _x1=0,_x2=overview_width - 1;
          int _y1,_y2;
          _y1=_y2=overlay_row[i] * overview_scale;
          XDrawLine(disp, pixmap, gc_select, _x1, _y1, _x2, _y2);
        }
        for(i = 0; i < cols; i++)
        {
          int _y1=0,_y2=overview_height - 1;
          int _x1,_x2;
          _x1 =_x2 =overlay_col[i] * overview_scale;
          XDrawLine(disp, pixmap, gc_select, _x1, _y1, _x2, _y2);
        }


/*
          for(i = 0; i < rows; i++)
          {
            XDrawLine(disp, pixmap, gc_select, 0, (overlay_row[i] - 1) * 
                      overview_scale, 
                      overview_width - 1, (overlay_row[i] - 1) * overview_scale);
          }
          for(i = 0; i < cols; i++)
          {
            XDrawLine(disp, pixmap, gc_select, (overlay_col[i] - 1) * 
                      overview_scale, 0, 
                      (overlay_col[i] - 1) * overview_scale, overview_height - 1);
          }
*/
          XSetForeground(disp, gc_select, SELECT_COLOR);

          XCopyArea(disp, pixmap, win, gc_select, 0, 0, overview_width, 
                  overview_height, 0, 0); 

          left_edge = (int)(GA_MIN(x1, x2) / overview_scale + .5);
          right_edge = (int)(GA_MAX(x1, x2) / overview_scale + .5);
          top_edge = (int)(GA_MIN(y1, y2) / overview_scale + .5);
          bottom_edge = (int)(GA_MAX(y1, y2) / overview_scale + .5);

          if (left_edge < 0) left_edge = 0;
          if (right_edge > grid_x - 1) right_edge = grid_x - 1;
          if (top_edge < 0) top_edge = 0;
          if (bottom_edge > grid_y - 1) bottom_edge = grid_y - 1;

          pict_width = right_edge - left_edge + 1; 
          pict_height = bottom_edge - top_edge + 1; 
          scale = GA_MIN( ceil((500.0 / pict_width)), 
                       ceil((500.0 / pict_height)));

          if ((pict_width < 501) || (pict_height < 501))
          {
            fprintf(stderr, "pict_width, pict_height, scale: %d, %d, %d\n",
                    pict_width, pict_height, scale);
            fprintf(stderr, "Left, right, top, bottom: %d, %d, %d, %d\n",
                    left_edge, right_edge, top_edge, bottom_edge);
            XtVaSetValues(view_button, XtNsensitive, True, NULL);
          }
          else
          {
            /* select to large a region to view - Error */
            XtVaSetValues(view_button, XtNsensitive, False, NULL);

            sprintf(err_str, "Selection Error: select area 500x500 or smaller");
            dialog_box(err_str);

            XSetForeground(disp, gc_select, CANVAS_COLOR);
            XDrawRectangle(disp, pixmap, gc_select, x1, y1, x2 - x1, y2 - y1); 
            XSetForeground(disp, gc_select, SELECT_COLOR);
            
            button_state = 0; 
            break;
          }
        }
      }
      break;

    case MotionNotify:
      if (button_state == 1)
      {
        if (another_state)
        {
        XSetForeground(disp, gc_select, CANVAS_COLOR);
        XFillRectangle(disp, pixmap, gc_select, 0, 0,overview_width, overview_height ); 
        XSetForeground(disp, gc_select, GRID_COLOR);
/*
        for(i = 0; i < rows; i++)
        {
          int x1=0,x2=overview_width - 1;
          int y1=y2=overlay_row[i] * overview_scale;
          XDrawLine(disp, pixmap, gc_select, x1, y1, x2, y2);
        }
        for(i = 0; i < cols; i++)
        {
          int y1=0,y2=overview_height - 1;
          int x1=x2=overlay_col[i] * overview_scale;
          XDrawLine(disp, pixmap, gc_select, x1, y1, x2, y2);
        }
*/

        for(i = 0; i < rows; i++)
        {
          XDrawLine(disp, pixmap, gc_select, 0, (overlay_row[i] - 0) * 
                    overview_scale, 
                    overview_width - 1, (overlay_row[i] - 0) * overview_scale);
        }
        for(i = 0; i < cols; i++)
        {
          XDrawLine(disp, pixmap, gc_select, (overlay_col[i] - 0) * 
                    overview_scale, 0, 
                    (overlay_col[i] - 0) * overview_scale, overview_height - 1);
        }
        XSetForeground(disp, gc_select, SELECT_COLOR);
        XCopyArea(disp, pixmap, win, gc_select, 0, 0, overview_width, 
                  overview_height, 0, 0); 
        XFlush(disp);
          another_state = False;
        }

        /* erase old rubberband box */
        XDrawRectangle(disp, win, gc_rband, x1, y1, x2 - x1, y2 - y1); 

        /* get latest corner */
        x2 = event->xbutton.x;
        y2 = event->xbutton.y;

        if (x2 < 0) x2 = 0;
        if (x2 > overview_width - 1) x2 = overview_width - 1;
        if (y2 < 0) y2 = 0;
        if (y2 > overview_height - 1) y2 = overview_height - 1;

        /* draw new rubberband box */
        XDrawRectangle(disp, win, gc_rband, x1, y1, x2 - x1, y2 - y1); 
      }
      break;
  } 

}

