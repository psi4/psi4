#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "xregion.h"
#include <math.h>

#define GA_ABS(x) (((x) < 0 )? -(x) : (x))

void DisplayIntervalValue()
{
  (void) sprintf(interval_string, "%4d ms", interval);
  XtSetArg(arg[0], XtNlabel, interval_string);
  XtSetValues(interval_widget,arg,1);
}

void DisplaySlowdownValue()
{
  (void) sprintf(slowdown_string, "%5.1f times", slowdown);
  XtSetArg(arg[0], XtNlabel, slowdown_string);
  XtSetValues(slowdown_widget,arg,1);
}

/**/
/* JJU: void ScrollProc(Widget scrollbar, caddr_t data, caddr_t position) */
void ScrollProc(Widget scrollbar, XtPointer data, XtPointer position)
/*
  Called when the left or right buttons are used to step the
  scrollbar left or right. We have the responsibility of
  moving the scrollbar.
*/
{
  Dimension length;
  float fraction;
  float shown;

  /* Get the scrollbar length and move the scroll bar */

  XtSetArg(arg[0], XtNlength, &length);
  XtGetValues(scrollbar, arg, 1);
  fraction = ((int) position)/ (double) length;

  interval = fraction*0.05*interval_max;   
  interval = GA_MIN(interval, interval_max);
  interval = GA_MAX(interval, 1);

  fraction = (float) interval/ (float) interval_max;
  shown = -1.0;

  DisplayIntervalValue();
  XawScrollbarSetThumb(scrollbar, fraction, shown);
}

/***** slowdown **********/
/* JJU: void ScrollProc2(Widget scrollbar, caddr_t data, caddr_t position) */
void ScrollProc2(Widget scrollbar, XtPointer data, XtPointer position)
/*
  Called when the left or right buttons are used to step the
  scrollbar left or right. We have the responsibility of
  moving the scrollbar.
*/
{
  Dimension length;
  double fraction;
  float shown;

  /* Get the scrollbar length and move the scroll bar */

  XtSetArg(arg[0], XtNlength, &length);
  XtGetValues(scrollbar, arg, 1);

  fraction -= ((int) position)/ (double) length;

  slowdown = fraction*0.2*slowdown_max;   


  /* need to add small number to avoid domain error in log(0) */

  slowdown = GA_MIN(slowdown,slowdown_max);
  slowdown = GA_MAX(slowdown, slowdown_min);


  /* scale current time according to the slowdown factor */
#ifdef DEBUG
  printf("before scaling %lu ( %ld %ld factor=%f) ",cur_time,
          slowdown,oldslowdown,(1.0*slowdown)/oldslowdown);
#endif /* DEBUG */
  cur_time = cur_time*slowdown/oldslowdown;

#ifdef DEBUG
    printf("and after %lu\n ",cur_time);
#endif /* DEBUG */
  oldslowdown = slowdown;

  fraction = (float) slowdown/ (float) slowdown_max;

  shown = -1.0;

  DisplaySlowdownValue();
  XawScrollbarSetThumb(scrollbar, (float)fraction, shown);
}
 

/**/
/* JJU: void JumpProc(Widget scrollbar, caddr_t data, caddr_t fraction_ptr) */
void JumpProc(Widget scrollbar, XtPointer data, XtPointer fraction_ptr) 
/*
  Called when the middle button is used to drag to 
  the scrollbar. The scrollbar is moved for us.
*/
{
  float fraction = *(float *) fraction_ptr;

  interval = fraction*interval_max;
  interval = GA_MIN(interval, interval_max);
  interval = GA_MAX(interval, 1);

  DisplayIntervalValue();
}

/**** slowdown ****/
/* JJU: void JumpProc2(Widget scrollbar, caddr_t data, caddr_t fraction_ptr) */
void JumpProc2(Widget scrollbar, XtPointer data, XtPointer fraction_ptr)
/*
  Called when the middle button is used to drag to 
  the scrollbar. The scrollbar is moved for us.
*/
{
  double exp_fraction;
  float fraction = *(float *) fraction_ptr;


  exp_fraction = pow(POW_BASE,(double)fraction);

  exp_fraction = (exp_fraction-1.)/ (pow(POW_BASE,1.) -1.);

  slowdown = exp_fraction*(slowdown_max-slowdown_min) + slowdown_min;
  slowdown = GA_MIN(slowdown,slowdown_max);
  slowdown = GA_MAX(slowdown, slowdown_min);


  /* scale current time according to the slowdown factor */
#ifdef DEBUG
  printf("before scaling %lu ( %ld %ld factor=%f) ",cur_time,
          slowdown,oldslowdown,(1.0*slowdown)/oldslowdown);
#endif /* DEBUG */
  cur_time = cur_time*slowdown/oldslowdown;

#ifdef DEBUG
  printf("and after %lu\n ",cur_time);
#endif /* DEBUG */
  oldslowdown = slowdown;

  DisplaySlowdownValue();
}
