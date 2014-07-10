#include <sys/types.h>
#include <netinet/in.h>
#include <stdio.h>
#include <string.h>
 
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/StringDefs.h>
#include <X11/Intrinsic.h>
#include <X11/IntrinsicP.h>
#include <X11/Shell.h>
#include <X11/ShellP.h>
#include <X11/Xaw/Scrollbar.h>
#include <X11/Xaw/Form.h>
#include <X11/Xaw/Command.h>
#include <X11/Xaw/Label.h>

/* #define DEBUG */
#define INITSLOW 50.
#define POW_BASE 100.

extern void exit();

/* function prototypes */
extern void Error(char *message, int err_num);
extern void UpdatePixRegion(int ilo, int ihi, int jlo, int jhi, int increment, 
                    double time);
extern void DisplayPixRegion(int ilo, int ihi, int jlo, int jhi);
extern void DisplaySlowdownValue();
extern void DisplayIntervalValue();
/* JJU: extern void ScrollProc(Widget scrollbar, caddr_t data,
   caddr_t position); */
extern void ScrollProc(Widget scrollbar, XtPointer data, XtPointer position);
/* JJU: extern void ScrollProc2(Widget scrollbar, caddr_t data,
   caddr_t position); */
extern void ScrollProc2(Widget scrollbar, XtPointer data, XtPointer position);
/* JJU: extern void JumpProc(Widget scrollbar, caddr_t data,
   caddr_t fraction_ptr); */
extern void JumpProc(Widget scrollbar, XtPointer data, XtPointer fraction_ptr);
/* JJU: extern void JumpProc2(Widget scrollbar, caddr_t data,
   caddr_t fraction_ptr); */
extern void JumpProc2(Widget scrollbar, XtPointer data,
                      XtPointer fraction_ptr);
extern void DrawColorMap();
extern void PrintColorMapText();
extern void UpdateDisplay();
/* JJU: extern void TimeOutCallback(caddr_t data); */
extern void TimeOutCallback(XtPointer data, XtIntervalId *xtintervalid);
/* JJU: extern void Exposed(Widget widget, caddr_t data, XEvent *event); */
extern void Exposed(Widget widget, XtPointer data, XEvent *event,
                    Boolean *bln);
/* JJU: extern void dismiss_dialog(Widget w, caddr_t data, XEvent *event) */
extern void dismiss_dialog(Widget widget, XtPointer data, XtPointer event);
/* JJU: extern void Quit(Widget widget, caddr_t data, XEvent *event); */
extern void Quit(Widget widget, XtPointer data, XtPointer event);
/* JJU: extern void StartStop(Widget widget, caddr_t data, XEvent *event); */
extern void StartStop(Widget widget, XtPointer data, XtPointer event);
/* JJU: extern void start_view(Widget widget, caddr_t data, XEvent *event); */
extern void start_view(Widget widget, XtPointer data, XtPointer event);
extern void Setcmap();
extern void ReadEventFile(char *filename);
/* JJU: extern void running_coords(Widget widget, caddr_t data,
   XEvent *event); */
extern void running_coords(Widget widget, XtPointer data, XEvent *event,
                           Boolean *bln);
/* JJU: extern void running_overview(Widget widget, caddr_t data,
   XEvent *event); */
extern void running_overview(Widget widget, XtPointer data, XEvent *event,
                             Boolean *bln);
/* JJU: extern void draw_select_box(Widget widget, caddr_t data,
   XEvent *event); */
extern void draw_select_box(Widget widget, XtPointer data, XEvent *event,
                            Boolean *bln);
extern int in_display_region(int ilo, int ihi, int jlo, int jhi);
extern XtAppContext create_overview(int argc, char **argv);
extern void create_main_window();
extern void setup_drawing();
extern void dialog_box(char *dlg_str);
extern void set_config();

/* Globals needed for display etc. */
#ifdef FIRST_TIME 
#define SCOPE 
#else
#define SCOPE extern
#endif /* FIRST_TIME */

SCOPE  Widget top_level_widget, box_widget, start_stop_button,
               scroll_widget, interval_widget,
               scroll_widget2, slowdown_widget,
               coord_widget, dlg_top, dlg_btn,  
               dlg_form, dlg_label, view_button,
               interval_label, slowdown_label,
               map_widget, quit_button, 
               overview_widget, select_widget,
               overview_title, overview_shell,
               canvas_widget, title_widget;
SCOPE  XtAppContext xregion_app;
SCOPE  XtIntervalId timer;
SCOPE  long first_time; /* Used to set scroll bar on first expose */

SCOPE  long interval_max;
SCOPE  long interval;         /* 0.5s between exposures by default */

SCOPE  double slowdown_max;
SCOPE  double slowdown_min;
SCOPE  double slowdown;        /* slowdown factor for animation */
SCOPE  double oldslowdown;

SCOPE  unsigned long int cur_time; /* current time */

SCOPE  Arg arg[25]; 
SCOPE  Display *display;
SCOPE  Window window, window_map;
SCOPE  int screen, depth;
SCOPE  Visual *visual;
SCOPE  XImage *image;
SCOPE  u_char *pict;
SCOPE  GC gc, gc_map;
SCOPE  char title[80];
SCOPE  char interval_string[10], slowdown_string[11];
SCOPE  int top_edge, bottom_edge, left_edge, right_edge;

#define GA_MAX(a,b) (((a)>(b)) ? (a) : (b))
#define GA_MIN(a,b) (((a)<(b)) ? (a) : (b))

#define MAX_COL 16
SCOPE u_char cmap[MAX_COL+1];
SCOPE Colormap colormap;

SCOPE int grid_x, grid_y;	/* The size of the grid */
SCOPE int scale;	/* No. of pixels per element */
SCOPE int pict_width;	/* The size of the picture = grid_x * scale */
SCOPE int pict_height;	/* The size of the picture = grid_y * scale */
SCOPE int overview_height, overview_width;
SCOPE double overview_scale; 
SCOPE int *overlay_row, *overlay_col;
SCOPE int rows, cols;
SCOPE  u_char *grid;

SCOPE  double *ltime;          /* last event time */
SCOPE  double *integr;         /* access integral */
SCOPE  double maxval;       /* max value of integral, zero is default */
SCOPE  u_char *flag;           /* access flag */

SCOPE  int working, animation;

/*** trace variables and constants ***/

SCOPE  long int num_events;
SCOPE  int *record;			/* tracefile data */
SCOPE  unsigned long int *ev_times;	/* times of events */
SCOPE  int cur_event;
#define RECLEN 8

/*** end of trace variables and constants ***/

/* color constants */
SCOPE Pixel DEFAULT_FG;
SCOPE Pixel DEFAULT_BG;
SCOPE Pixel RBAND_COLOR;
SCOPE Pixel SELECT_COLOR;
SCOPE Pixel GRID_COLOR;
SCOPE Pixel CANVAS_COLOR;

