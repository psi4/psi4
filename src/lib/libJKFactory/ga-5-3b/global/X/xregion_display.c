#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "xregion.h"

void UpdateDisplay()
{
int n;
char string[60];
XFontStruct *labelfont;
char *fontname = "8x13bold";
Window W;

  /* Incorporate the display changes after the animation is over */

  /* change the title */

  strcpy(title, "Time Lost Due To Contention");
  XtVaSetValues(title_widget,XtNlabel, title, NULL);

  /* delete all irrelevant widgets */

  XtDestroyWidget(slowdown_widget);
  XtDestroyWidget(scroll_widget2);
  XtDestroyWidget(interval_widget); 
  XtDestroyWidget(scroll_widget);
  XtDestroyWidget(start_stop_button);
  XtDestroyWidget(interval_label);
  XtDestroyWidget(slowdown_label);
  
  /* load the label font */
  
  if ((labelfont = XLoadQueryFont(display,fontname)) == NULL)
  {
     Error("failed to load label font",-1);
  }
       
  /* Add label widget to display the max value of integrals */

  sprintf(string, "max value = %8.5f", maxval);
  interval_widget = XtVaCreateManagedWidget("maxval", labelWidgetClass,
                                           box_widget,
                                           XtNbackground, DEFAULT_BG,
                                           XtNforeground, DEFAULT_FG,
                                           XtNx, 25,
                                           XtNwidth, 300,
                                           XtNvertDistance, 5,
                                           XtNfromVert, title_widget,
                                           XtNlabel, string,
                                           XtNborderWidth, 0,
                                           NULL); 

  XtVaSetValues(interval_widget, XtNfont, labelfont, NULL);

  /* Raise Quit widget partially obscurred by centered max value string */

  W = XtWindow(quit_button);
  XRaiseWindow(display,W);
 
  /* Now, update the ColorMap legend */

  XClearWindow(display, window_map);
  DrawColorMap();
  PrintColorMapText();

  XFlush(display);
}


/*****************************************************************/
/* JJU: void TimeOutCallback(caddr_t data) */
void TimeOutCallback(XtPointer data, XtIntervalId *xtintervalid)
{
#define  AMP 1
int ilo, ihi, jlo, jhi, inc;
int xlo, xhi, ylo, yhi;   /* display coordinates */
int base, stime;

/* Do work on time out here */

/*
   printf("time %lu, event=%d\n",cur_time,cur_event);
   fflush(stdout);
   printf("ev_times[cur_event]= ");
*/
  while (slowdown * (ev_times[cur_event] / 1000) < (cur_time + interval) &&
         cur_event < num_events)
  {
    base = cur_event*RECLEN;
    ilo = record[base+2];
    ihi = record[base+3];
    jlo = record[base+4];
    jhi = record[base+5];
    inc = AMP*record[base+7];
 
    if (in_display_region(ilo, ihi, jlo, jhi))
    {
      ilo = GA_MAX(ilo - top_edge, 0);
      ihi = GA_MIN(ihi - top_edge, pict_height - 1);
      jlo = GA_MAX(jlo - left_edge, 0);
      jhi = GA_MIN(jhi - left_edge, pict_width - 1);
      UpdatePixRegion(ilo, ihi, jlo, jhi, inc, (double) 1e-6 * ev_times[cur_event]);
      DisplayPixRegion(ilo, ihi, jlo, jhi);
    } 
    cur_event++;
  } 

  if(cur_event < num_events)
  {
    cur_time += interval;
  }
  else if(animation)
  {
    animation = False;
    printf("\nEnd of Event Animation ...\n");
    UpdatePixRegion(0, pict_height - 1, 0, pict_width - 1,  0, 0.0);
    DisplayPixRegion( 0, pict_height - 1, 0, pict_width - 1);
    UpdateDisplay();

#ifdef DEBUG
    for(inc = 0; inc < pict_width * pict_height; inc++)
    {
      printf("%f\n", integr[inc]);
    }
    printf("Intialized or computed max integral value = %f\n",maxval); 
#endif /* DEBUG */
  }  

  /* Restore the call back  */
  timer = XtAppAddTimeOut(xregion_app, interval, TimeOutCallback, NULL);
}  


/**/
/* JJU: void Exposed(Widget widget, caddr_t data, XEvent *event) */
void Exposed(Widget widget, XtPointer data, XEvent *event, Boolean *bln)
{
    /* Now we are exposed so we can draw ... */
    
    if (event->xexpose.count == 0)
    {
        if (first_time) 
        {
            /* Cannot seem to set this before now ? */
            ScrollProc(scroll_widget, NULL, 0);
            ScrollProc2(scroll_widget2, NULL, 0);
            first_time = False;
        }
        
        DrawColorMap();
        PrintColorMapText();
        DisplayPixRegion(0, pict_height - 1, 0, pict_width - 1);
        XFlush(display);
    }
}  

